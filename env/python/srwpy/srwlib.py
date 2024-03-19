#############################################################################
# SRWLib for Python v 0.20
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility

try: #OC15112022
    from . import srwlpy as srwl
    from . import uti_math
    from .srwl_uti_cryst import *
    from .uti_math_eigen import UtiMathEigen
except: #OC15112022
    import srwlpy as srwl
    import uti_math
    from srwl_uti_cryst import *
    from uti_math_eigen import UtiMathEigen

#import srwlpy as srwl
from array import *
from math import *
from copy import *

import datetime
import json
import random
import sys
import os
import traceback
#import uti_math
import errno
import tempfile
import shutil
import time
#import gc #OC13122023 TEST

#from srwl_uti_cryst import *
#from uti_math_eigen import UtiMathEigen #OC21062021

#try:
#    from uti_plot import * #universal simple plotting module distributed together with SRWLib
#except:
#    #excInf = sys.exc_info()
#    #print(excInf[1]) #printing exception value
#    traceback.print_exc()
#    print('Plotting utilities module was not loaded.')
#    print('1D and 2D plotting (generation of graphs, image plots, etc.) will not be possible.')

#****************************************************************************
#****************************************************************************
# Global Constants
#****************************************************************************
#****************************************************************************

_Pi = 3.14159265358979
_ElCh = 1.60217646263E-19 #1.602189246E-19 #Electron Charge [Q]
_ElMass_kg = 9.1093818872E-31 #9.10953447E-31 #Electron Mass in [kg]
_ElMass_MeV = 0.51099890221 #Electron Mass in [MeV]
_LightSp = 2.9979245812E+08 #Speed of Light [m/c]
_Light_eV_mu = 1.23984197 #OC23062019 #1.23984186 #Wavelength <-> Photon Energy conversion constant ([um] <-> [eV])
_PlanckConst_eVs = 4.13566766225E-15 #Planck constant in [eV*s]

#****************************************************************************
#****************************************************************************
# SRWLib Python Classes
#****************************************************************************
#****************************************************************************
class SRWLParticle(object):
    """Charged Particle"""

    def __init__(self, _x=0, _y=0, _z=0, _xp=0, _yp=0, _gamma=1, _relE0=1, _nq=-1):
        """
        :param _x: instant coordinates [m]
        :param _y: instant coordinates [m]
        :param _z: instant coordinates [m]
        :param _xp: instant transverse velocity component btx = vx/c (angles for relativistic particle)
        :param _yp: instant transverse velocity component bty = vy/c (angles for relativistic particle)
        :param _gamma: relative energy
        :param _relE0: rest mass (energy) in units of electron rest mass, e.g. 1 for electron, 1836.1526988 (=938.272013/0.510998902) for proton
        :param _nq: charge of the particle related to absolute value of electron charge, -1 for electron, 1 for positron and for proton
        """
        self.x = _x
        self.y = _y
        self.z = _z
        self.xp = _xp
        self.yp = _yp
        self.gamma = _gamma
        self.relE0 = _relE0
        self.nq = _nq

    def drift(self, _dist):
        """Propagates particle beam statistical moments over a distance in free space
        :param _dist: distance the beam has to be propagated over [m]
        """
        self.z += _dist
        self.x +=  self.xp*_dist
        self.y +=  self.yp*_dist
        #print('e-beam positions after drift: x=', self.x, ' y=', self.y, ' at z=', self.z)

    def get_E(self, _unit='GeV'):
        en = self.gamma*self.relE0*_ElMass_MeV #[MeV]
        if _unit == 'TeV': en *= 1e-06
        elif _unit == 'GeV': en *= 1e-03
        elif _unit == 'keV': en *= 1e+03
        elif _unit == 'eV': en *= 1e+06
        elif _unit == 'meV': en *= 1e+09
        return en

#****************************************************************************
class SRWLPartBeam(object):
    """Particle Beam"""

    def __init__(self, _Iavg=0, _nPart=0, _partStatMom1=None, _arStatMom2=None):
        """
        :param _Iavg: average current [A]
        :param _nPart: number of electrons (in a bunch)
        :param _partStatMom1: particle type and 1st order statistical moments
        :param _arStatMom2: 2nd order statistical moments
            [0]: <(x-x0)^2>
            [1]: <(x-x0)*(xp-xp0)>
            [2]: <(xp-xp0)^2>
            [3]: <(y-y0)^2>
            [4]: <(y-y0)*(yp-yp0)>
            [5]: <(yp-yp0)^2>
            [6]: <(x-x0)*(y-y0)>
            [7]: <(xp-xp0)*(y-y0)>
            [8]: <(x-x0)*(yp-yp0)>
            [9]: <(xp-xp0)*(yp-yp0)>
            [10]: <(E-E0)^2>/E0^2
            [11]: <(s-s0)^2>
            [12]: <(s-s0)*(E-E0)>/E0
            [13]: <(x-x0)*(E-E0)>/E0
            [14]: <(xp-xp0)*(E-E0)>/E0
            [15]: <(y-y0)*(E-E0)>/E0
            [16]: <(yp-yp0)*(E-E0)>/E0
            [17]: <(x-x0)*(s-s0)>
            [18]: <(xp-xp0)*(s-s0)>
            [19]: <(y-y0)*(s-s0)>
            [20]: <(yp-yp0)*(s-s0)>
        """
        self.Iavg = _Iavg
        self.nPart = _nPart
        self.partStatMom1 = SRWLParticle() if _partStatMom1 is None else _partStatMom1
        self.arStatMom2 = array('d', [0] * 21) if _arStatMom2 is None else _arStatMom2

    def from_Twiss(self, _Iavg=0, _e=0, _sig_e=0, _emit_x=0, _beta_x=0, _alpha_x=0, _eta_x=0, _eta_x_pr=0, _emit_y=0, _beta_y=0, _alpha_y=0, _eta_y=0, _eta_y_pr=0):
        """Sets up particle (electron) beam internal data from Twiss parameters
        :param _Iavg: average current [A]
        :param _e: energy [GeV]
        :param _sig_e: RMS energy spread
        :param _emit_x: horizontal emittance [m]
        :param _beta_x: horizontal beta-function [m]
        :param _alpha_x: horizontal alpha-function [rad]
        :param _eta_x: horizontal dispersion function [m]
        :param _eta_x_pr: horizontal dispersion function derivative [rad]
        :param _emit_y: vertical emittance [m]
        :param _beta_y: vertical beta-function [m]
        :param _alpha_y: vertical alpha-function [rad]
        :param _eta_y: vertical dispersion function [m]
        :param _eta_y_pr: vertical dispersion function derivative [rad]
        """
        self.Iavg = _Iavg
        self.partStatMom1.gamma = _e/0.51099890221e-03 #assuming electrons
        sigeE2 = _sig_e*_sig_e
        self.arStatMom2[0] = _emit_x*_beta_x + sigeE2*_eta_x*_eta_x #<(x-<x>)^2>
        self.arStatMom2[1] = -_emit_x*_alpha_x + sigeE2*_eta_x*_eta_x_pr #<(x-<x>)(x'-<x'>)>
        self.arStatMom2[2] = _emit_x*(1 + _alpha_x*_alpha_x)/_beta_x + sigeE2*_eta_x_pr*_eta_x_pr #<(x'-<x'>)^2>
        self.arStatMom2[3] = _emit_y*_beta_y + sigeE2*_eta_y*_eta_y #<(y-<y>)^2>
        self.arStatMom2[4] = -_emit_y*_alpha_y + sigeE2*_eta_y*_eta_y_pr #<(y-<y>)(y'-<y'>)>
        self.arStatMom2[5] = _emit_y*(1 + _alpha_y*_alpha_y)/_beta_y + sigeE2*_eta_y_pr*_eta_y_pr #<(y'-<y'>)^2>
        self.arStatMom2[10] = sigeE2

    def from_RMS(self, _Iavg=0, _e=0, _sig_e=0, _sig_x=0, _sig_x_pr=0, _m_xx_pr=0, _sig_y=0, _sig_y_pr=0, _m_yy_pr=0):
        """Sets up particle (electron) beam internal data from Twiss parameters
        :param _Iavg: average current [A]
        :param _e: energy [GeV]
        :param _sig_e: RMS energy spread
        :param _sig_x: horizontal RMS size [m]
        :param _sig_x_pr: horizontal RMS divergence [rad]
        :param _m_xx_pr: <(x-<x>)(x'-<x'>)> [m]
        :param _sig_y: vertical RMS size [m]
        :param _sig_y_pr: vertical RMS divergence [rad]
        :param _m_yy_pr: <(y-<y>)(y'-<y'>)> [m]
        """
        self.Iavg = _Iavg
        self.partStatMom1.gamma = _e/0.51099890221e-03 #assuming electrons
        sigeE2 = _sig_e*_sig_e
        self.arStatMom2[0] = _sig_x*_sig_x #<(x-<x>)^2>
        self.arStatMom2[1] = _m_xx_pr #<(x-<x>)(x'-<x'>)>
        self.arStatMom2[2] = _sig_x_pr*_sig_x_pr #<(x'-<x'>)^2>
        self.arStatMom2[3] = _sig_y*_sig_y #<(y-<y>)^2>
        self.arStatMom2[4] = _m_yy_pr #<(y-<y>)(y'-<y'>)>
        self.arStatMom2[5] = _sig_y_pr*_sig_y_pr #<(y'-<y'>)^2>
        self.arStatMom2[10] = sigeE2
        
    def drift(self, _dist):
        """Propagates particle beam statistical moments over a distance in free space
        :param _dist: distance the beam has to be propagated over [m]
        """
        self.partStatMom1.drift(_dist)
        self.arStatMom2[0] += self.arStatMom2[1]*_dist*2 + self.arStatMom2[2]*_dist*_dist
        self.arStatMom2[1] += self.arStatMom2[2]*_dist
        self.arStatMom2[3] += self.arStatMom2[4]*_dist*2 + self.arStatMom2[5]*_dist*_dist
        self.arStatMom2[4] += self.arStatMom2[5]*_dist
        #print('E-beam 2nd order moments after drift=', _dist)
        #print(self.arStatMom2)
        #to be checked and extended for other stat. moments

#****************************************************************************
class SRWLMagFld(object):
    """Magnetic Field (base class)"""
    
class SRWLMagFld3D(SRWLMagFld):
    """Magnetic Field: Arbitrary 3D"""
    
    def __init__(self, _arBx=None, _arBy=None, _arBz=None, _nx=0, _ny=0, _nz=0, _rx=0, _ry=0, _rz=0, _nRep=1, _interp=1, _arX=None, _arY=None, _arZ=None):
        """
        :param _arBx: horizontal magnetic field component array [T]
        :param _arBy: vertical magnetic field component array [T]
        :param _arBz: longitudinal magnetic field component array [T]
        :param _nx: number of magnetic field data points in the horizontal direction
        :param _ny: number of magnetic field data points in the vertical direction
        :param _nz: number of magnetic field data points in the longitudinal direction
        :param _rx: range of horizontal coordinate for which the field is defined [m]
        :param _ry: range of vertical coordinate for which the field is defined [m]
        :param _rz: range of longitudinal coordinate for which the field is defined [m]
        :param _nRep: "number of periods", i.e. number of times the field is "repeated" in the longitudinal direction
        :param _interp: interpolation method to use (e.g. for trajectory calculation), 1- bi-linear (3D), 2- (bi-)quadratic (3D), 3- (bi-)cubic (3D)
        :param _arX: optional array of horizontal transverse coordinate of an irregular 3D mesh (if this array is defined, rx will be ignored)
        :param _arY: optional array of vertical transverse coordinate of an irregular 3D mesh (if this array is defined, ry will be ignored)
        :param _arZ: optional array of longitudinal coordinate of an irregular 3D mesh (if this array is defined, rz will be ignored)
        """
        self.arBx = array('d') if _arBx is None else _arBx
        self.arBy = array('d') if _arBy is None else _arBy
        self.arBz = array('d') if _arBz is None else _arBz
        self.nx = _nx
        self.ny = _ny
        self.nz = _nz
        self.rx = _rx
        self.ry = _ry
        self.rz = _rz
        self.arX = array('d') if _arX is None else _arX
        self.arY = array('d') if _arY is None else _arY
        self.arZ = array('d') if _arZ is None else _arZ
        self.nRep = _nRep
        self.interp = _interp

    def add_const(self, _bx=0, _by=0, _bz=0):
        """Adds constant magnetic field to the entire tabulated field (to simulate e.g. background magnetic field effects)
        :param _bx: horizontal magnetic field component to add [T]
        :param _by: vertical magnetic field component to add [T]
        :param _bz: longitudinal magnetic field component to add [T]
        """
        nTot = self.nx*self.ny*self.nz
        for i in range(nTot):
            self.arBx[i] += _bx
            self.arBy[i] += _by
            self.arBz[i] += _bz

    def save_ascii(self, _file_path, _xc=0, _yc=0, _zc=0):
        """Auxiliary function to write tabulated Arbitrary 3D Magnetic Field data to ASCII file"""
        sHead = '#Bx [T], By [T], Bz [T] on 3D mesh: inmost loop vs X (horizontal transverse position), outmost loop vs Z (longitudinal position)\n'
        sHead += '#' + repr(-0.5*self.rx + _xc) + ' #initial X position [m]\n'
        sHead += '#' + repr(0. if(self.nx <= 1) else self.rx/(self.nx - 1)) + ' #step of X [m]\n'
        sHead += '#' + repr(self.nx) + ' #number of points vs X\n'
        sHead += '#' + repr(-0.5*self.ry + _yc) + ' #initial Y position [m]\n'
        sHead += '#' + repr(0. if(self.ny <= 1) else self.ry/(self.ny - 1)) + ' #step of Y [m]\n'
        sHead += '#' + repr(self.ny) + ' #number of points vs Y\n'
        sHead += '#' + repr(-0.5*self.rz + _zc) + ' #initial Z position [m]\n'
        sHead += '#' + repr(0. if(self.nz <= 1) else self.rz/(self.nz - 1)) + ' #step of Z [m]\n'
        sHead += '#' + repr(self.nz) + ' #number of points vs Z\n'
        arColsWr = [self.arBx, self.arBy, self.arBz]
        #print(self.nx, self.rx, self.ny, self.ry, self.nz, self.rz)
        srwl_uti_write_data_cols(_file_path, arColsWr, '\t', sHead)

class SRWLMagFldM(SRWLMagFld):
    """Magnetic Field: Multipole Magnet"""
    
    #def __init__(self, _G=0, _m=2, _n_or_s='n', _Leff=0, _Ledge=0):
    def __init__(self, _G=0, _m=2, _n_or_s='n', _Leff=0, _Ledge=0, _R=0):
        """
        :param _G: field parameter [T] for dipole, [T/m] for quadrupole (negative means defocusing for x), [T/m^2] for sextupole, [T/m^3] for octupole
        :param _m: multipole order 1 for dipole, 2 for quadrupoole, 3 for sextupole, 4 for octupole
        :param _n_or_s: normal ('n') or skew ('s')
        :param _Leff: effective length [m]
        :param _Ledge: "soft" edge length for field variation from 10% to 90% [m]; G/(1 + ((z-zc)/d)^2)^2 fringe field dependence is assumed
        :param _R: radius of curvature of central trajectory [m] (for simulating e.g. quadrupole component integrated to a bending magnet; effective if > 0)
        """
        self.G = _G
        self.m = _m
        self.n_or_s = _n_or_s
        self.Leff = _Leff
        self.Ledge = _Ledge
        self.R = _R

class SRWLMagFldS(SRWLMagFld):
    """Magnetic Field: Solenoid"""
    
    def __init__(self, _B=0, _Leff=0):
        """
        :param _B: magnetic field [T]
        :param _Leff: effective length [m]
        """
        self.B = _B
        self.Leff = _Leff

class SRWLMagFldH(SRWLMagFld):
    """Magnetic Field: Undulator Harmonic"""
    
    def __init__(self, _n=1, _h_or_v='v', _B=0, _ph=0, _s=1, _a=1):
        """
        :param _n: harmonic number
        :param _h_or_v: magnetic field plane horzontal ('h') or vertical ('v')
        :param _B: magnetic field amplitude [T]
        :param _ph: initial phase [rad]
        :param _s: symmetry vs longitudinal position 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
        :param _a: coefficient for transverse depenednce B*cosh(2*Pi*n*a*y/per)*cos(2*Pi*n*z/per + ph)
        """
        self.n = _n
        self.h_or_v = _h_or_v
        self.B = _B
        self.ph = _ph
        self.s = _s
        self.a = _a

class SRWLMagFldU(SRWLMagFld):
    """Magnetic Field: Undulator"""
    
    def __init__(self, _arHarm=None, _per=0, _nPer=0):
        """
        :param _arHarm: array of field harmonics
        :param _per: period length [m]
        :param _nPer: number of periods (will be rounded to integer)
        """
        self.arHarm = [] if _arHarm is None else _arHarm
        self.per = _per
        self.nPer = _nPer

    def allocate(self, _nHarm):
        #self.arHarm = [SRWLMagFldH()]*_nHarm
        #arHarmLoc = []
        #for i in range(_nHarm): arHarm.append(SRWLMagFldH())
        self.arHarm = [] #OC18112019
        for i in range(_nHarm): self.arHarm.append(SRWLMagFldH())

    def set_sin(self, _per=0.02, _len=1, _bx=0, _by=0, _phx=0, _phy=0, _sx=1, _sy=1):
        """Setup basic undulator with sinusoidal magnetic field
        :param _per: period length [m]
        :param _len: undulator length [m]
        :param _bx: horizontal magnetic field amplitude [m]
        :param _by: vertical magnetic field amplitude [m]
        :param _phx: initial phase of the horizontal magnetic field [rad]
        :param _phx: initial phase of the vertical magnetic field [rad]
        :param _sx: symmetry of the horizontal magnetic field vs longitudinal position 1 - symmetric (B ~ cos(2*Pi*n*z/per + ph)) , -1 - anti-symmetric (B ~ sin(2*Pi*n*z/per + ph))
        :param _sy: symmetry of the vertical magnetic field vs longitudinal position
        """
        nPerAvg = int(round(_len/_per))
        self.nPer = nPerAvg
        self.per = _per
        if(len(self.arHarm) > 0):
            del self.arHarm; self.arHarm = []
        if(_bx != 0): self.arHarm.append(SRWLMagFldH(_h_or_v='h', _B=_bx, _ph=_phx, _s=_sx))
        if(_by != 0): self.arHarm.append(SRWLMagFldH(_h_or_v='v', _B=_by, _ph=_phy, _s=_sy))

    def get_K(self):
        """Estimate K (deflection parameter) value"""
        mult = _ElCh/(2.*_Pi*_ElMass_kg*_LightSp)
        nHarm = len(self.arHarm)
        sumBdNe2 = 0
        for i in range(nHarm):
            curHarm = self.arHarm[i]
            curBdN = curHarm.B/curHarm.n
            sumBdNe2 += curBdN*curBdN
        return mult*self.per*sqrt(sumBdNe2)

    def K_2_B(self, K): #MR31072016 (added)
        """Convert K (deflection parameter) to B (magnetic field amplitude)"""
        mult = _ElCh/(2.*_Pi*_ElMass_kg*_LightSp)
        B = K / (mult * self.per)
        return B

    def get_E1(self, _en_elec=3., _unit='eV'):
        """Estimate fundamental photon energy
        :param _en_elec: electron energy [GeV]
        :return: fundamental photon energy [eV]
        """
        K = self.get_K()
        gamma =  1000.*_en_elec/_ElMass_MeV
        lamda_m = self.per*(1. + 0.5*K*K)/(2.*gamma*gamma)
        return srwl_uti_ph_en_conv(lamda_m, _in_u='m', _out_u=_unit)

    def E1_2_K(self, _e1, _en_elec=3.):
        """Estimate deflection parameter from 
        :param _e1: fundamental photon energy [eV]
        :param _en_elec: electron energy [GeV]
        :return: deflection parameter
        """
        buf = 9.4963421866853*_en_elec*_en_elec/self.per/_e1
        if(buf < 1): return 0
        else: return sqrt((buf - 1)*2)

    def E1_2_B(self, _e1, _en_elec=3.):
        """Estimate deflection parameter from 
        :param _e1: fundamental photon energy [eV]
        :param _en_elec: electron energy [GeV]
        :return: magnetic field amplitude [T]
        """
        K = self.E1_2_K(_e1, _en_elec)
        return 2*_Pi*_ElMass_kg*_LightSp*K/(_ElCh*self.per)

class SRWLMagFldC(SRWLMagFld):
    """Magnetic Field: Container"""
    
    def __init__(self, _arMagFld=None, _arXc=None, _arYc=None, _arZc=None, _arVx=None, _arVy=None, _arVz=None, _arAng=None):
    #def __init__(self, _arMagFld=None, _arXc=None, _arYc=None, _arZc=None):
        """
        :param _arMagFld: magnetic field structures array
        :param _arXc: horizontal center positions of magnetic field elements in arMagFld array [m]
        :param _arYc: vertical center positions of magnetic field elements in arMagFld array [m]
        :param _arZc: longitudinal center positions of magnetic field elements in arMagFld array [m]
        :param _arVx: horizontal components of axis vectors of magnetic field elements in arMagFld array [rad]
        :param _arVy: vertical components of axis vectors of magnetic field elements in arMagFld array [rad]
        :param _arVz: longitudinal components of axis vectors of magnetic field elements in arMagFld array [rad]
        :param _arAng: rotation angles of magnetic field elements about their axes [rad]
        """
        #self.arMagFld = [] if _arMagFld is None else _arMagFld

        if(_arMagFld is None):
            self.arMagFld = []
            self.arXc = array('d') if _arXc is None else _arXc
            self.arYc = array('d') if _arYc is None else _arYc
            self.arZc = array('d') if _arZc is None else _arZc
            #The following arrays are optional
            #self.arVx = array('d') if _arVx is None else _arVx
            #self.arVy = array('d') if _arVy is None else _arVy
            #self.arVz = array('d') if _arVz is None else _arVz
            #self.arAng = array('d') if _arAng is None else _arAng
        else:
            if(not(isinstance(_arMagFld, list) or isinstance(_arMagFld, array) or isinstance(_arMagFld, tuple))):
                self.arMagFld = [_arMagFld] #to allow for simple initialization by one element
                nElem = 1
            else:
                self.arMagFld = _arMagFld
                nElem = len(_arMagFld)

            if(_arXc is None):
                self.arXc = array('d', [0]*nElem)
            elif(isinstance(_arXc, array)):
                self.arXc = _arXc
            elif(isinstance(_arXc, list)): #or isinstance(_arXc, tuple)):
                self.arXc = array('d', _arXc)
            elif(nElem == 1):
                self.arXc = array('d', [0])
                self.arXc[0] = _arXc

            if(_arYc is None):
                self.arYc = array('d', [0]*nElem)
            #elif(isinstance(_arYc, list) or isinstance(_arYc, array) or isinstance(_arYc, tuple)):
            #    self.arYc = _arYc
            elif(isinstance(_arYc, array)):
                self.arYc = _arYc
            elif(isinstance(_arYc, list)): #or isinstance(_arYc, tuple)):
                self.arYc = array('d', _arYc)
            elif(nElem == 1):
                self.arYc = array('d', [0])
                self.arYc[0] = _arYc

            if(_arZc == None):
                self.arZc = array('d', [0]*nElem)
            #elif(isinstance(_arZc, list) or isinstance(_arZc, array) or isinstance(_arZc, tuple)):
            #    self.arZc = _arZc
            elif(isinstance(_arZc, array)):
                self.arZc = _arZc
            elif(isinstance(_arZc, list)): #or isinstance(_arZc, tuple)):
                self.arZc = array('d', _arZc)
            elif(nElem == 1):
                self.arZc = array('d', [0])
                self.arZc[0] = _arZc

            arVxWasSubm = False
            if(_arVx is None):
                self.arVx = array('d', [0]*nElem)
            elif(isinstance(_arVx, array)):
                self.arVx = _arVx
                arVxWasSubm = True
            elif(isinstance(_arVx, list)):
                self.arVx = array('d', _arVx)
                arVxWasSubm = True
            elif(nElem == 1):
                self.arVx = array('d', [0])
                self.arVx[0] = _arVx

            arVyWasSubm = False
            if(_arVy is None):
                self.arVy = array('d', [0]*nElem)
            elif(isinstance(_arVy, array)):
                self.arVy = _arVy
                arVyWasSubm = True
            elif(isinstance(_arVy, list)):
                self.arVy = array('d', _arVy)
                arVyWasSubm = True
            elif(nElem == 1):
                self.arVy = array('d', [0])
                self.arVy[0] = _arVy
                
            if(_arVz is None):
                self.arVz = array('d', [1]*nElem)
                if(arVxWasSubm and arVyWasSubm):
                    lenArVx = len(_arVx)
                    lenArVy = len(_arVy)
                    if(lenArVx == lenArVy):
                        for i in range(lenArVx):
                            self.arVz[i] = sqrt(1. - _arVx[i]*_arVx[i] - _arVy[i]*_arVy[i])
            elif(isinstance(_arVz, array)):
                self.arVz = _arVz
            elif(isinstance(_arVz, list)):
                self.arVz = array('d', _arVz)
            elif(nElem == 1):
                self.arVz = array('d', [1])
                self.arVz[0] = _arVz

            if(_arAng is None):
                self.arAng = array('d', [0]*nElem)
            elif(isinstance(_arAng, array)):
                self.arAng = _arAng
            elif(isinstance(_arAng, list)):
                self.arAng = array('d', _arAng)
            elif(nElem == 1):
                self.arAng = array('d', [0])
                self.arAng[0] = _arAng

    def allocate(self, _nElem):
        self.arMagFld = [SRWLMagFld()]*_nElem
        self.arXc = array('d', [0]*_nElem)
        self.arYc = array('d', [0]*_nElem)
        self.arZc = array('d', [0]*_nElem)
        self.arVx = array('d', [0]*_nElem)
        self.arVy = array('d', [0]*_nElem)
        self.arVz = array('d', [1]*_nElem)
        self.arAng = array('d', [0]*_nElem)

    def add(self, _mag, _xc=None, _yc=None, _zc=None, _vx=None, _vy=None, _vz=None, _ang=None):
        """Adds magnetic element to container
        :param _mag: magnetic element (or array of elements) to be added
        :param _xc: horizontal center position (or array of center positions) of magnetic field element to be added [m]
        :param _yc: vertical center positions (or array of center positions) of magnetic field element to be added [m]
        :param _zc: longitudinal center positions (or array of center positions) of magnetic field element to be added [m]
        :param _vx: horizontal component of axis vectors of magnetic field element to be added [rad]
        :param _vy: vertical component of axis vectors of magnetic field element to be added [rad]
        :param _vz: longitudinal components of axis vector of magnetic field element to be added [rad]
        :param _ang: rotation angle about axis [rad]
        """
        if(_mag is None):
            raise Exception("No magnetic field elements were supplied for adding to container") 
        if(isinstance(_mag, list) or isinstance(_mag, array)):
            lenMag = len(_mag)
            if((_xc is None) and (_yc is None) and (_zc is None) and
               (_vx is None) and (_vy is None) and (_vz is None) and (_ang is None)):
                for i in range(lenMag): self.add(_mag[i])
            elif((isinstance(_xc, list) or isinstance(_xc, array)) and
                 (isinstance(_yc, list) or isinstance(_yc, array)) and
                 (isinstance(_zc, list) or isinstance(_zc, array)) and
                 (isinstance(_vx, list) or isinstance(_vx, array)) and
                 (isinstance(_vy, list) or isinstance(_vy, array)) and
                 (isinstance(_vz, list) or isinstance(_vz, array)) and
                 (isinstance(_ang, list) or isinstance(_ang, array))):
                lenXc = len(_xc)
                lenYc = len(_yc)
                lenZc = len(_zc)
                lenVx = len(_vx)
                lenVy = len(_vy)
                lenVz = len(_vz)
                lenAng = len(_ang)
                if((lenXc == lenMag) and (lenYc == lenMag) and (lenZc == lenMag) and
                   (lenVx == lenMag) and (lenVy == lenMag) and (lenVz == lenMag) and (lenAng == lenMag)):
                    for i in range(lenMag): self.add(_mag[i], _xc[i], _yc[i], _zc[i], _vx[i], _vy[i], _vz[i])
                else: raise Exception("Inconsistent magnetic element positions data") 
        else:
            self.arMagFld.append(_mag)
            if(_xc is None): _xc = 0
            if(_yc is None): _yc = 0
            if(_zc is None): _zc = 0
            if(_vx is None): _vx = 0
            if(_vy is None): _vy = 0
            if(_vz is None):
                _vz = 1.
                if((_vx is not None) and (_vy is not None)):
                    _vz = sqrt(1. - _vx*_vx - _vy*_vy)
            if(_ang is None): _ang = 0
            self.arXc.append(_xc)
            self.arYc.append(_yc)
            self.arZc.append(_zc)
            self.arVx.append(_vx)
            self.arVy.append(_vy)
            self.arVz.append(_vz)
            self.arAng.append(_ang)

#****************************************************************************
class SRWLPrtTrj(object):
    """Charged Particle Trajectory"""

    def __init__(self, _arX=None, _arXp=None, _arY=None, _arYp=None, _arZ=None, _arZp=None, _arBx=None, _arBy=None, _arBz=None, _np=0, _ctStart=0, _ctEnd=0, _partInitCond=None):
        """
        :param _arX: array of horizontal position [m]
        :param _arXp: array of horizontal relative velocity (trajectory angle) [rad]
        :param _arY: array of vertical position [m]
        :param _arYp: array of vertical relative velocity (trajectory angle) [rad]
        :param _arZ: array of longitudinal positions [m]
        :param _arZp: array of longitudinal relative velocity [rad]
        :param _arBx: array of horizontal magnetic field component "seen" by particle [T]
        :param _arBy: array of vertical magnetic field component "seen" by particle [T]
        :param _arBz: array of longitudinal magnetic field component "seen" by particle [T]
        :param _np: number of trajectory points
        :param _ctStart: start value of independent variable (c*t) for which the trajectory should be (/is) calculated (is constant step enough?)
        :param _ctEnd: end value of independent variable (c*t) for which the trajectory should be (/is) calculated (is constant step enough?)
        :param _partInitCond: particle type and initial conditions for which the trajectory should be (/is) calculated
        """

        if(_np > 0):
            self.arX = array('d', [0]*_np) if _arX is None else _arX
            self.arY = array('d', [0]*_np) if _arY is None else _arY
            self.arZ = array('d', [0]*_np) if _arZ is None else _arZ
            self.arXp = array('d', [0]*_np) if _arXp is None else _arXp
            self.arYp = array('d', [0]*_np) if _arYp is None else _arYp
            self.arZp = array('d', [0]*_np) if _arZp is None else _arZp
        else:
            self.arX = array('d') if _arX is None else _arX
            self.arY = array('d') if _arY is None else _arY
            self.arZ = array('d') if _arZ is None else _arZ
            self.arXp = array('d') if _arXp is None else _arXp
            self.arYp = array('d') if _arYp is None else _arYp
            self.arZp = array('d') if _arZp is None else _arZp
            
        if _arBx is not None: self.arBx = _arBx #by default, arBx, _arBy, arBz are not created
        if _arBy is not None: self.arBy = _arBy
        if _arBz is not None: self.arBz = _arBz
            
        self.np = _np
        self.ctStart = _ctStart
        self.ctEnd = _ctEnd
        self.partInitCond = SRWLParticle() if _partInitCond is None else _partInitCond

    def allocate(self, _np, _allB=False):
        _np = int(_np)
        self.arX = array('d', [0]*_np)
        self.arXp = array('d', [0]*_np)
        self.arY = array('d', [0]*_np)
        self.arYp = array('d', [0]*_np)
        self.arZ = array('d', [0]*_np)
        self.arZp = array('d', [0]*_np)
        self.np = _np
        if _allB == True:
            self.arBx = array('d', [0]*_np)
            self.arBy = array('d', [0]*_np)
            self.arBz = array('d', [0]*_np)

    def save_ascii(self, _file_path):
        """Auxiliary function to write tabulated Trajectory data to ASCII file"""
        f = open(_file_path, 'w')
        resStr = '#ct [m], X [m], BetaX [rad], Y [m], BetaY [rad], Z [m], BetaZ [rad]'
        if(hasattr(self, 'arBx')):
            resStr += ', Bx [T]'
        if(hasattr(self, 'arBy')):
            resStr += ', By [T]'
        if(hasattr(self, 'arBz')):
            resStr += ', Bz [T]'
        f.write(resStr + '\n')
        ctStep = 0
        if self.np > 0:
            ctStep = (self.ctEnd - self.ctStart)/(self.np - 1)
        ct = self.ctStart
        for i in range(self.np):
            resStr = str(ct) + '\t' + repr(self.arX[i]) + '\t' + repr(self.arXp[i]) + '\t' + repr(self.arY[i]) + '\t' + repr(self.arYp[i]) + '\t' + repr(self.arZ[i]) + '\t' + repr(self.arZp[i])
            if(hasattr(self, 'arBx')):
                resStr += '\t' + repr(self.arBx[i])
            if(hasattr(self, 'arBy')):
                resStr += '\t' + repr(self.arBy[i])
            if(hasattr(self, 'arBz')):
                resStr += '\t' + repr(self.arBz[i])
            f.write(resStr + '\n')        
            ct += ctStep
        f.close()
      
#****************************************************************************
class SRWLKickM(object):
    """Kick Matrix (for fast trajectory calculation)"""
    
    def __init__(self, _arKickMx=None, _arKickMy=None, _order=2, _nx=0, _ny=0, _nz=0, _rx=0, _ry=0, _rz=0, _x=0, _y=0, _z=0):
        """
        :param _arKickMx: horizontal kick-matrix (tabulated on the same transverse grid vs x and y as vertical kick-matrix)
        :param _arKickMy: vertical kick-matrix (tabulated on the same transverse grid vs x and y as horizontal kick-matrix)
        :param _order: kick order: 1- first order (in this case kick matrix data is assumed to be in [T*m]), 2- second order (kick matrix data is assumed to be in [T^2*m^2])
        :param _nx: numbers of points in kick matrices in horizontal direction
        :param _ny: numbers of points in kick matrices in vertical direction
        :param _nz: number of steps in longitudinal direction
        :param _rx: range covered by kick matrices in horizontal direction [m]
        :param _ry: range covered by kick matrices in vertical direction [m]
        :param _rz: extension in longitudinal direction [m]
        :param _x: horizontal coordinate of center point [m]
        :param _y: vertical coordinate of center point [m]
        :param _z: longitudinal coordinate of center point [m]
        """
        self.arKickMx = array('d') if _arKickMx is None else _arKickMx
        self.arKickMy = array('d') if _arKickMy is None else _arKickMy
        self.order = _order
        self.nx = _nx
        self.ny = _ny
        self.nz = _nz
        self.rx = _rx
        self.ry = _ry
        self.rz = _rz
        self.x = _x
        self.y = _y
        self.z = _z

#****************************************************************************
class SRWLGsnBm(object):
    """(Coherent) Gaussian (Radiation) Beam"""
    
    def __init__(self, _x=0, _y=0, _z=0, _xp=0, _yp=0, _avgPhotEn=1, _pulseEn=1, _repRate=1, _polar=1, _sigX=10e-06,
                 _sigY=10e-06, _sigT=1e-15, _mx=0, _my=0):
                 #_sigY=10e-06, _sigT=1e-15, _mx=0, _my=0, _arIncohStatMom2=None): #OC15092017 (?)
        """
        :param _x: average horizontal coordinate of waist [m]
        :param _y: average vertical coordinate of waist [m]
        :param _z: average longitudinal coordinate of waist [m]
        :param _xp: average horizontal angle at waist [rad]
        :param _yp: average verical angle at waist [rad]
        :param _avgPhotEn: average photon energy [eV]
        :param _pulseEn: energy per pulse [J]
        :param _repRate: rep. rate [Hz]
        :param _polar: polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left
        :param _sigX: rms beam size vs horizontal position [m] at waist (for intensity)
        :param _sigY: rms beam size vs vertical position [m] at waist (for intensity)
        :param _sigT: rms pulse duration [s] (for intensity)
        :param _mx: transverse Gauss-Hermite mode order in horizontal direction
        :param _my: transverse Gauss-Hermite mode order in vertical direction
        """
        self.x = _x
        self.y = _y
        self.z = _z
        self.xp = _xp
        self.yp = _yp
        self.avgPhotEn = _avgPhotEn
        self.pulseEn = _pulseEn
        self.repRate = _repRate
        self.polar = _polar
        self.sigX = _sigX
        self.sigY = _sigY
        self.sigT = _sigT
        self.mx = _mx
        self.my = _my

        #if(_arIncohStatMom2 is not None): self.arIncohStatMom2 = _arIncohStatMom2 #OC15092017 (?)

#****************************************************************************
class SRWLPtSrc(object):
    """Point Source (emitting coherent spherical wave)"""
    
    def __init__(self, _x=0, _y=0, _z=0, _flux=1, _unitFlux=1, _polar=1):
        """
        :param _x: horizontal position [m]
        :param _y: vertical position [m]
        :param _z: longitudinal position [m]
        :param _flux: spectral flux value
        :param _unitFlux: spectral flux units: 1- ph/s/.1%bw, 2- W/eV
        :param _polar: polarization 1- lin. hor., 2- lin. vert., 3- lin. 45 deg., 4- lin.135 deg., 5- circ. right, 6- circ. left, 7- radial
        """
        self.x = _x
        self.y = _y
        self.z = _z
        self.flux = _flux
        self.unitFlux = _unitFlux
        self.polar = _polar

#****************************************************************************
class SRWLRadMesh(object):
    """Radiation Mesh (Sampling)"""
    
    def __init__(self, _eStart=0, _eFin=0, _ne=1, _xStart=0, _xFin=0, _nx=1, _yStart=0, _yFin=0, _ny=1, _zStart=0, _nvx=0, _nvy=0, _nvz=1, _hvx=1, _hvy=0, _hvz=0, _arSurf=None):
        """
        :param _eStart: initial value of photon energy (/time)
        :param _eFin: final value of photon energy (/time)
        :param _ne: number of points vs photon energy (/time)
        :param _xStart: initial value of horizontal position (/angle)
        :param _xFin: final value of horizontal position (/angle)
        :param _nx: number of points vs horizontal position (/angle)
        :param _yStart: initial value of vertical position (/angle)
        :param _yFin: final value of vertical position (/angle)
        :param _ny: number of points vs vertical position (/angle)
        :param _zStart: longitudinal position
        :param _nvx: horizontal lab-frame coordinate of inner normal to observation plane (/ surface in its center)
        :param _nvy: vertical lab-frame coordinate of inner normal to observation plane (/ surface in its center)
        :param _nvz: longitudinal lab-frame coordinate of inner normal to observation plane (/ surface in its center)
        :param _hvx: horizontal lab-frame coordinate of the horizontal base vector of the observation plane (/ surface in its center)
        :param _hvy: vertical lab-frame coordinate of the horizontal base vector of the observation plane (/ surface in its center)
        :param _hvz: longitudinal lab-frame coordinate of the horizontal base vector of the observation plane (/ surface in its center)
        :param _arSurf: array defining the observation surface (as function of 2 variables - x & y - on the mesh given by _xStart, _xFin, _nx, _yStart, _yFin, _ny; to be used in case this surface differs from plane)
        """
        self.eStart = _eStart
        self.eFin = _eFin
        self.ne = _ne
        self.xStart = _xStart
        self.xFin = _xFin
        self.nx = _nx
        self.yStart = _yStart
        self.yFin = _yFin
        self.ny = _ny
        self.zStart = _zStart

        self.nvx = _nvx
        self.nvy = _nvy
        self.nvz = _nvz
        self.hvx = _hvx
        self.hvy = _hvy
        self.hvz = _hvz
        self.arSurf = _arSurf

    def set_from_other(self, _mesh):
        self.eStart = _mesh.eStart; self.eFin = _mesh.eFin; self.ne = _mesh.ne;
        self.xStart = _mesh.xStart; self.xFin = _mesh.xFin; self.nx = _mesh.nx;
        self.yStart = _mesh.yStart; self.yFin = _mesh.yFin; self.ny = _mesh.ny;
        self.zStart = _mesh.zStart

        self.nvx = _mesh.nvx; self.nvy = _mesh.nvy; self.nvz = _mesh.nvz
        self.hvx = _mesh.hvx; self.hvy = _mesh.hvy; self.hvz = _mesh.hvz
        del self.arSurf; self.arSurf = None
        if(_mesh.arSurf is not None):
            try:
                lenArSurf = len(_mesh.arSurf)
                if(lenArSurf > 0):
                    self.arSurf = array('d', [0]*lenArSurf)
                    for i in range(lenArSurf):
                        self.arSurf[i] = _mesh.arSurf[i]
            except:
                pass

    def is_equal(self, _mesh, _rel_tol=1.e-07, _check_surf=False):
        #isEqual = True
        if(_mesh.ne != self.ne): return False
        if(_mesh.nx != self.nx): return False
        if(_mesh.ny != self.ny): return False

        absTolE = _rel_tol*self.eStart
        if((self.eFin > 0) and (self.eFin != self.eStart) and (self.ne > 1)): absTolE = _rel_tol*abs((self.eFin - self.eStart)/(self.ne - 1))
        if(abs(self.eStart - _mesh.eStart) > absTolE): return False
        if((self.eFin > 0) and (_mesh.eFin > 0)):
            if(abs(self.eFin - _mesh.eFin) > absTolE): return False
        
        absTolX = _rel_tol*self.xStart
        if((self.xFin > 0) and (self.xFin != self.xStart) and (self.nx > 1)): absTolX = _rel_tol*abs((self.xFin - self.xStart)/(self.nx - 1))
        if(absTolX == 0): absTolX = 1.e-13 #to tune
        if(abs(self.xStart - _mesh.xStart) > absTolX): return False
        if((self.xFin > 0) and (_mesh.xFin > 0)):
            if(abs(self.xFin - _mesh.xFin) > absTolX): return False

        absTolY = _rel_tol*self.yStart
        if((self.yFin > 0) and (self.yFin != self.yStart) and (self.ny > 1)): absTolY = _rel_tol*abs((self.yFin - self.yStart)/(self.ny - 1))
        if(absTolY == 0): absTolY = 1.e-13 #to tune
        if(abs(self.yStart - _mesh.yStart) > absTolY): return False
        if((self.yFin > 0) and (_mesh.yFin > 0)):
            if(abs(self.yFin - _mesh.yFin) > absTolY): return False

        if(abs(_mesh.nvx - self.nvx) > _rel_tol): return False
        if(abs(_mesh.nvy - self.nvy) > _rel_tol): return False
        if(abs(_mesh.nvz - self.nvz) > _rel_tol): return False
        if(abs(_mesh.hvx - self.hvx) > _rel_tol): return False
        if(abs(_mesh.hvy - self.hvy) > _rel_tol): return False
        if(abs(_mesh.hvz - self.hvz) > _rel_tol): return False

        if(not _check_surf): return True
        if(self.arSurf is None):
            if(_mesh.arSurf is None): return True
            else: return False
        else:
            if(_mesh.arSurf is None): return False
            else:
                lenArSurf = len(self.arSurf)
                if(lenArSurf != len(_mesh.arSurf)): return False
                for i in range(lenArSurf):
                    if(self.arSurf[i] != _mesh.arSurf[i]): return False

        return True

    def get_dep_type(self): #Get possible dependency type (for intensity calc.)
        depType = -1
        if((self.ne >= 1) and (self.nx == 1) and (self.ny == 1)): depType = 0
        elif((self.ne == 1) and (self.nx > 1) and (self.ny == 1)): depType = 1
        elif((self.ne == 1) and (self.nx == 1) and (self.ny > 1)): depType = 2
        elif((self.ne == 1) and (self.nx > 1) and (self.ny > 1)): depType = 3
        elif((self.ne > 1) and (self.nx > 1) and (self.ny == 1)): depType = 4
        elif((self.ne > 1) and (self.nx == 1) and (self.ny > 1)): depType = 5
        elif((self.ne > 1) and (self.nx > 1) and (self.ny > 1)): depType = 6
        return depType

    def copy(self): #is called from C++
        return deepcopy(self)

#****************************************************************************
class SRWLStokes(object):
    """Radiation Stokes Parameters"""
    
    #def __init__(self, _arS0=None, _arS1=None, _arS2=None, _arS3=None, _typeStokes='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0):
    #def __init__(self, _arS=None, _typeStokes='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0, _mutual=0):
    #def __init__(self, _arS=None, _typeStokes='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0, _mutual=0, _n_comp=4): #OC04022021
    def __init__(self, _arS=None, _typeStokes='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0, _mutual=0, _n_comp=4, _itStFin=None): #OC04022021
        """
        :param _arS: flat C-aligned array of all Stokes components (outmost loop over Stokes parameter number); NOTE: only 'f' (float) is supported for the moment (Jan. 2012)
        :param _typeStokes: numerical type: 'f' (float) or 'd' (double, not supported yet)
        :param _eStart: initial value of photon energy (/time)
        :param _eFin: final value of photon energy (/time)
        :param _ne: numbers of points vs photon energy
        :param _xStart: initial value of horizontal position
        :param _xFin: final value of photon horizontal position
        :param _nx: numbers of points vs horizontal position
        :param _yStart: initial value of vertical position
        :param _yFin: final value of vertical position
        :param _ny: numbers of points vs vertical position
        :param _mutual: switch specifying that mutual Stokes components should be or are defined (2*_n_comp*(_ne*_nx*_ny_)^2 values)
        :param _n_comp: number of Stoke components (1 to 4)
        :param _itStFin: pair of start and end index values of general conjugated coordinate (for partial allocation of Mutual Intensity)
        """
        self.arS = _arS #flat C-aligned array of all Stokes components (outmost loop over Stokes parameter number); NOTE: only 'f' (float) is supported for the moment (Jan. 2012)
        self.numTypeStokes = _typeStokes #electric field numerical type: 'f' (float) or 'd' (double)
        self.mesh = SRWLRadMesh(_eStart, _eFin, _ne, _xStart, _xFin, _nx, _yStart, _yFin, _ny) #to make mesh an instance variable
        self.avgPhotEn = 0 #average photon energy for time-domain simulations    
        self.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
        self.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
        self.unitStokes = 1 #Stokes units: 0- arbitrary, 1- Phot/s/0.1%bw/mm^2 ?
        self.mutual = _mutual #indicator of Mutual Stokes components

        nProd = _ne*_nx*_ny #array length to store one component of complex electric field
        if(nProd > 0): #OC29062021
            if(_arS is not None):
                if(isinstance(_arS, int)):
                    if(_arS == 1): self.allocate(_ne, _nx, _ny, _typeStokes, _mutual, _n_comp, _itStFin) #OC29062021
            
        #if((_arS == 1) and (nProd > 0)):
            #self.allocate(_ne, _nx, _ny, _typeStokes, _mutual, _n_comp, _itStFin) #OC03032021     
            #self.allocate(_ne, _nx, _ny, _typeStokes, _mutual, _n_comp) #OC04022021     
            #self.allocate(_ne, _nx, _ny, _typeStokes, _mutual)          
        #s0needed = 0
        #s1needed = 0
        #s2needed = 0
        #s3needed = 0
        #if((_arS0 == 1) and (nProd > 0)):
        #    s0needed = 1
        #if((_arS1 == 1) and (nProd > 0)):
        #    s1needed = 1
        #if((_arS2 == 1) and (nProd > 0)):
        #    s2needed = 1
        #if((_arS3 == 1) and (nProd > 0)):
        #    s3needed = 1
        #if((s0needed > 0) or (s1needed > 0) or (s2needed > 0) or (s3needed > 0)):
        #    self.allocate(_ne, _nx, _ny, s0needed, s1needed, s2needed, s3needed)

    #def allocate(self, _ne, _nx, _ny, s0needed=1, s1needed=1, s2needed=1, s3needed=1, _typeStokes='f'):
    #def allocate(self, _ne, _nx, _ny, _typeStokes='f', _mutual=0):
    #def allocate(self, _ne, _nx, _ny, _typeStokes='f', _mutual=0, _n_comp=4): #OC04022021
    def allocate(self, _ne, _nx, _ny, _typeStokes='f', _mutual=0, _n_comp=4, _itStFin=None): #OC03032021
        #print('') #debugging
        #print('          (re-)allocating: old point numbers: ne=',self.mesh.ne,' nx=',self.mesh.nx,' ny=',self.mesh.ny,' type:',self.numTypeStokes)
        #print('                           new point numbers: ne=',_ne,' nx=',_nx,' ny=',_ny,' type:',_typeStokes)
        #nTot = _ne*_nx*_ny #array length to one Stokes component
        #if s0needed:
        #    del self.arS0
        #    self.arS0 = array(_typeStokes, [0]*nTot)
        #if s1needed:
        #    del self.arS1
        #    self.arS1 = array(_typeStokes, [0]*nTot)
        #if s2needed:
        #    del self.arS2
        #    self.arS2 = array(_typeStokes, [0]*nTot)
        #if s3needed:
        #    del self.arS3
        #    self.arS3 = array(_typeStokes, [0]*nTot)

        nTot = _ne*_nx*_ny
        if _mutual > 0:
            #nTot *= nTot
            nTotIt = nTot #OC03032021
            if(_itStFin is not None): #OC03032021
                nTotIt = _ne*(_itStFin[1] - _itStFin[0] + 1)

            nTot *= (nTotIt*2) #OC03032021
            #nTot *= (nTot*2) #OC04052018 (since <E(r1)E*(r2)> may be a complex entity)
        #nTot *= 4 #array length of all Stokes components
        nTot *= _n_comp #OC04022021 #array length of all Stokes components
        #eventually allow for storage of less than 4 Stokes components!

        #DEBUG
        #print('_ne=', _ne, '_nx=', _nx, '_ny=', _ny, 'nTot=', nTot)
        #END DEBUG
        
        #self.arS = array(_typeStokes, [0]*nTot)
        #OC26022021: The above "allocation", passing through list, strongly "overconsumes" memory :(
        self.arS = srwl_uti_array_alloc(_typeStokes, nTot) #OC26022021
        #The above thing works much better (faster, less "overconsumption" of memory)
        #Do Python developers have a direct method of allocating arrays of given type without passing through lists?

        self.numTypeStokes = _typeStokes
        self.mesh.ne = _ne
        self.mesh.nx = _nx
        self.mesh.ny = _ny
        self.mutual = _mutual
        
        if(_mutual != 0): self.mesh.type = 2 #OC04032021 (this indicates MI in Python)

        if(_itStFin is not None): #OC03032021 (if this is in mesh, it can be easily submitted to a number of functions implemented in C++)
            self.mesh.itStart = _itStFin[0]
            self.mesh.itFin = _itStFin[1] 

    def add_stokes(self, _st, _n_comp=4, _mult=1, _meth=0):
        """Add Another Stokes structure
        :param _st: Stokes structure to be added
        :param _n_comp: number of components to treat
        :param _mult: multiplier 
        :param _meth: method of adding the Stokes structure _st:
        0- simple addition assuming _wfr to have same mesh as this wavefront
        1- add using bilinear interpolation (taking into account meshes of the two wavefronts)
        2- add using bi-quadratic interpolation (taking into account meshes of the two wavefronts)
        3- add using bi-cubic interpolation (taking into account meshes of the two wavefronts)
        """

        #nTot = _n_comp*self.mesh.ne*self.mesh.nx*self.mesh.ny #eventually allow for storage of less than 4 Stokes components!
        nTot = self.mesh.ne*self.mesh.nx*self.mesh.ny 
        if(self.mutual > 0): 
            #nTot *= nTot
            nTot *= (nTot*2) #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)

        nTot *= _n_comp #eventually allow for storage of less than 4 Stokes components!

        if(_meth == 0):
            if((self.mesh.ne != _st.mesh.ne) or (self.mesh.nx != _st.mesh.nx) or (self.mesh.ny != _st.mesh.ny)):
                raise Exception("Stokes parameters addition can not be performed by this method because of unequal sizes of the two Stokes structures") 

            st_arS = _st.arS
            if(_mult == 1):
                for i in range(nTot):
                    #for some reason, this increases memory requirements in Py(?):
                    self.arS[i] += st_arS[i]
            else:
                for i in range(nTot):
                    #for some reason, this increases memory requirements in Py(?):
                    self.arS[i] += _mult*st_arS[i]
                
        elif(_meth == 1):
            #to implement
            raise Exception("This Stokes parameters addition method is not implemented yet")
        
        elif(_meth == 2):
            #to implement
            raise Exception("This Stokes parameters addition method is not implemented yet")
            
        elif(_meth == 3):
            #to implement
            raise Exception("This Stokes parameters addition method is not implemented yet")

    def avg_update_same_mesh(self, _more_stokes, _iter, _n_stokes_comp=4, _mult=1., _sum=False): #OC20112020
    #def avg_update_same_mesh(self, _more_stokes, _iter, _n_stokes_comp=4, _mult=1.):
        """ Update this Stokes data structure with new data, contained in the _more_stokes structure, calculated on the same mesh, so that this structure would represent estimation of average of (_iter + 1) structures
        :param _more_stokes: Stokes data structure to "add" to the estimation of average
        :param _iter: number of Stokes structures already "added" previously
        :param _n_stokes_comp: number of Stokes components to treat (1 to 4)
        :param _mult: optional multiplier of the _more_stokes
        """

        #DEBUG
        #print('avg_update_same_mesh: iter=', _iter, _mult)

        #nStPt = self.mesh.ne*self.mesh.nx*self.mesh.ny*_n_stokes_comp
        nStPt = self.mesh.ne*self.mesh.nx*self.mesh.ny
        if self.mutual > 0:
            #nStPt *= nStPt
            nStPt *= (nStPt*2) #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)

        nStPt *= _n_stokes_comp

        if(_sum):
            if(_mult == 1.):
                for ir in range(nStPt): self.arS[ir] += _more_stokes.arS[ir]
            else:
                for ir in range(nStPt): self.arS[ir] += _mult*_more_stokes.arS[ir]
        else:
            if(_mult == 1.):
                for ir in range(nStPt):
                    self.arS[ir] = (self.arS[ir]*_iter + _more_stokes.arS[ir])/(_iter + 1)
            else:
                for ir in range(nStPt):
                    self.arS[ir] = (self.arS[ir]*_iter + _mult*_more_stokes.arS[ir])/(_iter + 1)

    def avg_update_interp(self, _more_stokes, _iter, _ord, _n_stokes_comp=4, _mult=1., _sum=False): #OC04112020
    #def avg_update_interp(self, _more_stokes, _iter, _ord, _n_stokes_comp=4, _mult=1.):
        """ Update this Stokes data structure with new data, contained in the _more_stokes structure, calculated on a different 2D mesh, so that it would represent estimation of average of (_iter + 1) structures
        :param _more_stokes: Stokes data structure to "add" to the estimation of average
        :param _iter: number of Stokes structures already "added" previously
        :param _ord: order of 2D interpolation to use (1- bilinear, ..., 3- bi-cubic)
        :param _n_stokes_comp: number of Stokes components to treat (1 to 4)
        :param _mult: optional multiplier of the _more_stokes
        """

        #DEBUG
        #print('avg_update_interp: iter=', _iter, _mult)
        #print('self.mesh.xStart=', self.mesh.xStart, 'self.mesh.xFin=', self.mesh.xFin, 'self.mesh.yStart=', self.mesh.yStart, 'self.mesh.yFin=', self.mesh.yFin)
        #print('_more_stokes.mesh.xStart=', _more_stokes.mesh.xStart, '_more_stokes.mesh.xFin=', _more_stokes.mesh.xFin, '_more_stokes.mesh.yStart=', _more_stokes.mesh.yStart, '_more_stokes.mesh.yFin=', _more_stokes.mesh.yFin)
        #END DEBUG
        
        eNpMeshRes = self.mesh.ne
        xNpMeshRes = self.mesh.nx
        xStartMeshRes = self.mesh.xStart
        xStepMeshRes = 0
        if(xNpMeshRes > 1):
            xStepMeshRes = (self.mesh.xFin - xStartMeshRes)/(xNpMeshRes - 1)
        yNpMeshRes = self.mesh.ny
        yStartMeshRes = self.mesh.yStart
        yStepMeshRes = 0
        if(yNpMeshRes > 1):
            yStepMeshRes = (self.mesh.yFin - yStartMeshRes)/(yNpMeshRes - 1)

        eNpWfr = _more_stokes.mesh.ne
        xStartWfr = _more_stokes.mesh.xStart
        xNpWfr = _more_stokes.mesh.nx
        #DEBUG
        #print('xNpWfr=', xNpWfr)
        #print('_more_stokes.mesh.xFin=', _more_stokes.mesh.xFin)
        #print('xStartWfr=', xStartWfr)

        xStepWfr = 0
        if(xNpWfr > 1):
            xStepWfr = (_more_stokes.mesh.xFin - xStartWfr)/(xNpWfr - 1)
        yStartWfr = _more_stokes.mesh.yStart
        yNpWfr = _more_stokes.mesh.ny
        yStepWfr = 0
        if(yNpWfr > 1):
            yStepWfr = (_more_stokes.mesh.yFin - yStartWfr)/(yNpWfr - 1)
        #DEBUG
        #print('avg_update_interp: iter=', _iter)
        #print('xStepWfr = ', xStepWfr)
        #END DEBUG

        nRadWfr = eNpWfr*xNpWfr*yNpWfr
        iOfstSt = 0
        ir = 0
        for iSt in range(_n_stokes_comp):
            for iy in range(yNpMeshRes):
                yMeshRes = yStartMeshRes + iy*yStepMeshRes
                for ix in range(xNpMeshRes):
                    xMeshRes = xStartMeshRes + ix*xStepMeshRes
                    for ie in range(eNpMeshRes):
                        #calculate Stokes parameters of propagated wavefront on the resulting mesh
                        #fInterp = srwl_uti_interp_2d(xMeshRes, yMeshRes, xStartWfr, xStepWfr, xNpWfr, yStartWfr, yStepWfr, yNpWfr, workArStokes, 1, eNpWfr, iOfstStokes)
                        fInterp = 0
                        loc_ix_ofst = iOfstSt + ie
                        nx_ix_per = xNpWfr*eNpWfr
                        if(_ord == 1): #bi-linear interpolation based on 4 points
                            ix0 = int(trunc((xMeshRes - xStartWfr)/xStepWfr + 1.e-09))
                            if((ix0 < 0) or (ix0 >= xNpWfr - 1)):
                                self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                continue
                            ix1 = ix0 + 1
                            tx = (xMeshRes - (xStartWfr + xStepWfr*ix0))/xStepWfr
                            iy0 = int(trunc((yMeshRes - yStartWfr)/yStepWfr + 1.e-09))
                            if((iy0 < 0) or (iy0 >= yNpWfr - 1)):
                                self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                continue
                            iy1 = iy0 + 1
                            ty = (yMeshRes - (yStartWfr + yStepWfr*iy0))/yStepWfr
                            iy0_nx_ix_per = iy0*nx_ix_per
                            iy1_nx_ix_per = iy1*nx_ix_per
                            ix0_ix_per_p_ix_ofst = ix0*eNpWfr + loc_ix_ofst
                            ix1_ix_per_p_ix_ofst = ix1*eNpWfr + loc_ix_ofst
                            a00 = _more_stokes.arS[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f10 = _more_stokes.arS[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            f01 = _more_stokes.arS[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f11 = _more_stokes.arS[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            a10 = f10 - a00
                            a01 = f01 - a00
                            a11 = a00 - f01 - f10 + f11
                            fInterp = a00 + tx*(a10 + ty*a11) + ty*a01

                            #DEBUG
                            #if(isnan(fInterp)):
                            #    print('nx=', _more_stokes.mesh.nx, ' ny=', _more_stokes.mesh.ny)
                            #    print('i00=', iy0_nx_ix_per + ix0_ix_per_p_ix_ofst, ' i10=', iy0_nx_ix_per + ix1_ix_per_p_ix_ofst, ' i01=', iy1_nx_ix_per + ix0_ix_per_p_ix_ofst, ' i11=', iy1_nx_ix_per + ix1_ix_per_p_ix_ofst)
                            #    print('ix0=', ix0, ' ix1=', ix1, ' iy0=', iy0, ' iy1=', iy1)
                            #    print('f00=', a00, ' f10=', f10, ' f01=', f01, ' f11=', f11)
                            #    sys.exit()
                            #else: print('OK')
                            #END DEBUG

                        elif(_ord == 2): #bi-quadratic interpolation based on 6 points
                            ix0 = int(round((xMeshRes - xStartWfr)/xStepWfr))
                            if((ix0 < 0) or (ix0 >= xNpWfr - 1)):
                                self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                continue
                            ixm1 = ix0 - 1
                            ix1 = ix0 + 1
                            tx = (xMeshRes - (xStartWfr + xStepWfr*ix0))/xStepWfr
                            iy0 = int(round((yMeshRes - yStartWfr)/yStepWfr))
                            if((iy0 < 0) or (iy0 >= yNpWfr - 1)):
                                self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                continue
                            iym1 = iy0 - 1
                            iy1 = iy0 + 1
                            ty = (yMeshRes - (yStartWfr + yStepWfr*iy0))/yStepWfr
                            iym1_nx_ix_per = iym1*nx_ix_per
                            iy0_nx_ix_per = iy0*nx_ix_per
                            iy1_nx_ix_per = iy1*nx_ix_per
                            ixm1_ix_per_p_ix_ofst = ixm1*eNpWfr + loc_ix_ofst
                            ix0_ix_per_p_ix_ofst = ix0*eNpWfr + loc_ix_ofst
                            ix1_ix_per_p_ix_ofst = ix1*eNpWfr + loc_ix_ofst
                            fm10 = _more_stokes.arS[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
                            a00 = _more_stokes.arS[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f10 = _more_stokes.arS[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            f0m1 = _more_stokes.arS[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f01 = _more_stokes.arS[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f11 = _more_stokes.arS[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            a10 = 0.5*(f10 - fm10)
                            a01 = 0.5*(f01 - f0m1)
                            a11 = a00 - f01 - f10 + f11
                            a20 = 0.5*(f10 + fm10) - a00
                            a02 = 0.5*(f01 + f0m1) - a00
                            fInterp = a00 + tx*(a10 + tx*a20 + ty*a11) + ty*(a01 + ty*a02)
    
                        elif(_ord == 3): #bi-cubic interpolation based on 12 points
                            ix0 = int(trunc((xMeshRes - xStartWfr)/xStepWfr + 1.e-09))
                            if((ix0 < 0) or (ix0 >= xNpWfr - 1)):
                                self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                continue
                            elif(ix0 < 1):
                                ix0 = 1
                            elif(ix0 >= xNpWfr - 2):
                                ix0 = xNpWfr - 3
                            ixm1 = ix0 - 1
                            ix1 = ix0 + 1
                            ix2 = ix0 + 2
                            tx = (xMeshRes - (xStartWfr + xStepWfr*ix0))/xStepWfr
                            iy0 = int(trunc((yMeshRes - yStartWfr)/yStepWfr + 1.e-09))
                            if((iy0 < 0) or (iy0 >= yNpWfr - 1)):
                                self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                continue
                            elif(iy0 < 1):
                                iy0 = 1
                            elif(iy0 >= yNpWfr - 2):
                                iy0 = yNpWfr - 3
                            iym1 = iy0 - 1
                            iy1 = iy0 + 1
                            iy2 = iy0 + 2
                            ty = (yMeshRes - (yStartWfr + yStepWfr*iy0))/yStepWfr
                            iym1_nx_ix_per = iym1*nx_ix_per
                            iy0_nx_ix_per = iy0*nx_ix_per
                            iy1_nx_ix_per = iy1*nx_ix_per
                            iy2_nx_ix_per = iy2*nx_ix_per
                            ixm1_ix_per_p_ix_ofst = ixm1*eNpWfr + loc_ix_ofst
                            ix0_ix_per_p_ix_ofst = ix0*eNpWfr + loc_ix_ofst
                            ix1_ix_per_p_ix_ofst = ix1*eNpWfr + loc_ix_ofst
                            ix2_ix_per_p_ix_ofst = ix2*eNpWfr + loc_ix_ofst
                            f0m1 = _more_stokes.arS[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f1m1 = _more_stokes.arS[iym1_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            fm10 = _more_stokes.arS[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
                            a00 = _more_stokes.arS[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f10 = _more_stokes.arS[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            f20 = _more_stokes.arS[iy0_nx_ix_per + ix2_ix_per_p_ix_ofst]
                            fm11 = _more_stokes.arS[iy1_nx_ix_per + ixm1_ix_per_p_ix_ofst]
                            f01 = _more_stokes.arS[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f11 = _more_stokes.arS[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            f21 = _more_stokes.arS[iy1_nx_ix_per + ix2_ix_per_p_ix_ofst]
                            f02 = _more_stokes.arS[iy2_nx_ix_per + ix0_ix_per_p_ix_ofst]
                            f12 = _more_stokes.arS[iy2_nx_ix_per + ix1_ix_per_p_ix_ofst]
                            a10 = -0.5*a00 + f10 - f20/6 - fm10/3
                            a01 = -0.5*a00 + f01 - f02/6 - f0m1/3
                            a11 = -0.5*(f01 + f10) + (f02 - f12 + f20 - f21)/6 + (f0m1 - f1m1 + fm10 - fm11)/3 + f11
                            a20 = -a00 + 0.5*(f10 + fm10)
                            a02 = -a00 + 0.5*(f01 + f0m1)
                            a21 = a00 - f01 + 0.5*(f11 - f10 - fm10 + fm11)
                            a12 = a00 - f10 + 0.5*(f11 - f01 - f0m1 + f1m1)
                            a30 = 0.5*(a00 - f10) + (f20 - fm10)/6
                            a03 = 0.5*(a00 - f01) + (f02 - f0m1)/6
                            a31 = 0.5*(f01 + f10 - f11 - a00) + (f21 + fm10 - f20 - fm11)/6
                            a13 = 0.5*(f10 - f11 - a00 + f01) + (f0m1 + f12 - f02 - f1m1)/6
                            fInterp = a00 + tx*(a10 + tx*(a20 + tx*(a30 + ty*a31) + ty*a21) + ty*a11) + ty*(a01 + ty*(a02 + ty*(a03 + tx*a13) + tx*a12))

                        #self.arS[ir] = (self.arS[ir]*_iter + fInterp)/(_iter + 1)
                        #self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)

                        #OC04112020
                        if(_sum): self.arS[ir] += _mult*fInterp
                        else: self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)
                        
                        ir += 1
            iOfstSt += nRadWfr        

    def avg_update_interp_mutual(self, _more_stokes, _iter, _n_stokes_comp=4, _mult=1., _sum=False): #OC13112020
    #def avg_update_interp_mutual(self, _more_stokes, _iter, _n_stokes_comp=4, _mult=1.):
        """ Update this Stokes data structure with new data, contained in the _more_stokes structure, calculated on a different 2D mesh, so that it would represent estimation of average of (_iter + 1) structures
        :param _more_stokes: Stokes data structure to "add" to the estimation of average
        :param _iter: number of Stokes structures already "added" previously
        :param _n_stokes_comp: number of Stokes components to treat (1 to 4)
        :param _mult: optional multiplier of the _more_stokes
        """
        eNpMeshRes = self.mesh.ne
        eStartMeshRes = self.mesh.eStart
        eStepMeshRes = 0
        if(eNpMeshRes > 1):
            eStepMeshRes = (self.mesh.eFin - eStartMeshRes)/(eNpMeshRes - 1)
        
        xNpMeshRes = self.mesh.nx
        xStartMeshRes = self.mesh.xStart
        xStepMeshRes = 0
        if(xNpMeshRes > 1):
            xStepMeshRes = (self.mesh.xFin - xStartMeshRes)/(xNpMeshRes - 1)
            
        yNpMeshRes = self.mesh.ny
        yStartMeshRes = self.mesh.yStart
        yStepMeshRes = 0
        if(yNpMeshRes > 1):
            yStepMeshRes = (self.mesh.yFin - yStartMeshRes)/(yNpMeshRes - 1)

        eNpWfr = _more_stokes.mesh.ne
        eStartWfr = _more_stokes.mesh.eStart
        eNpWfr = _more_stokes.mesh.ne
        eStepWfr = 0
        if(eNpWfr > 1):
            eStepWfr = (_more_stokes.mesh.eFin - eStartWfr)/(eNpWfr - 1)
        eNpWfr_mi_1 = eNpWfr - 1
        
        xStartWfr = _more_stokes.mesh.xStart
        xNpWfr = _more_stokes.mesh.nx
        xStepWfr = 0
        if(xNpWfr > 1):
            xStepWfr = (_more_stokes.mesh.xFin - xStartWfr)/(xNpWfr - 1)
        xNpWfr_mi_1 = xNpWfr - 1

        yStartWfr = _more_stokes.mesh.yStart
        yNpWfr = _more_stokes.mesh.ny
        yStepWfr = 0
        if(yNpWfr > 1):
            yStepWfr = (_more_stokes.mesh.yFin - yStartWfr)/(yNpWfr - 1)
        yNpWfr_mi_1 = yNpWfr - 1

        #nRadWfr = eNpWfr*xNpWfr*yNpWfr
        #nRadWfr *= nRadWfr

        #DEBUG
        #print('Base Mesh X:', self.mesh.xStart, 0 if(self.mesh.nx <= 1) else (self.mesh.xFin - self.mesh.xStart)/(self.mesh.nx - 1), self.mesh.nx)
        #print('More Mesh X:', xStartWfr, xStepWfr, xNpWfr)
        #print('Base Mesh Y:', self.mesh.yStart, 0 if(self.mesh.ny <= 1) else (self.mesh.yFin - self.mesh.yStart)/(self.mesh.ny - 1), self.mesh.ny)
        #print('More Mesh Y:', yStartWfr, yStepWfr, yNpWfr)

        perEp = 2 #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)
        perE = perEp*eNpWfr #OC04052018 
        #perE = eNpWfr
        perXp = perE*eNpWfr
        perX = perXp*xNpWfr
        perYp = perX*xNpWfr
        perY = perYp*yNpWfr
        nRadWfr = perY*yNpWfr
        
        iter_p_1 = _iter + 1 #OC04052018
        iter_d_iter_p_1 = _iter/iter_p_1

        iOfstSt = 0
        ir = 0
        for iSt in range(_n_stokes_comp):
            for iy in range(yNpMeshRes):
                doZeroFy = False
                yMeshRes = yStartMeshRes + iy*yStepMeshRes
                iy0 = 0
                if(yStepWfr > 0): iy0 = int(trunc((yMeshRes - yStartWfr)/yStepWfr + 1.e-09))
                if((iy0 < 0) or (iy0 > yNpWfr_mi_1)):
                    doZeroFy = True
                    #self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                    #continue
                iy1 = iy0 + 1
                if(iy1 > yNpWfr_mi_1): iy1 = yNpWfr_mi_1
                ty = 0
                if(yStepWfr > 0): ty = (yMeshRes - (yStartWfr + yStepWfr*iy0))/yStepWfr

                iy0_perY = iy0*perY
                iy1_perY = iy1*perY
                
                for iyp in range(yNpMeshRes):
                    doZeroFyp = False
                    ypMeshRes = yStartMeshRes + iyp*yStepMeshRes
                    iyp0 = 0
                    if(yStepWfr > 0): iyp0 = int(trunc((ypMeshRes - yStartWfr)/yStepWfr + 1.e-09))
                    if((iyp0 < 0) or (iyp0 > yNpWfr_mi_1)):
                        doZeroFyp = True
                        #self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                        #continue
                    iyp1 = iyp0 + 1
                    if(iyp1 > yNpWfr_mi_1): iyp1 = yNpWfr_mi_1
                    typ = 0
                    if(yStepWfr > 0): typ = (ypMeshRes - (yStartWfr + yStepWfr*iyp0))/yStepWfr

                    iyp0_perYp = iyp0*perYp
                    iyp1_perYp = iyp1*perYp
                    iyp0_perYp_p_iy0_perY = iyp0_perYp + iy0_perY
                    iyp1_perYp_p_iy0_perY = iyp1_perYp + iy0_perY
                    iyp0_perYp_p_iy1_perY = iyp0_perYp + iy1_perY
                    iyp1_perYp_p_iy1_perY = iyp1_perYp + iy1_perY

                    for ix in range(xNpMeshRes):
                        doZeroFx = False
                        xMeshRes = xStartMeshRes + ix*xStepMeshRes
                        ix0 = 0
                        if(xStepWfr > 0): ix0 = int(trunc((xMeshRes - xStartWfr)/xStepWfr + 1.e-09))
                        if((ix0 < 0) or (ix0 > xNpWfr_mi_1)):
                            doZeroFx = True
                            #self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                            #continue
                        ix1 = ix0 + 1
                        if(ix1 > xNpWfr_mi_1): ix1 = xNpWfr_mi_1
                        tx = 0
                        if(xStepWfr > 0): tx = (xMeshRes - (xStartWfr + xStepWfr*ix0))/xStepWfr

                        ix0_perX = ix0*perX
                        ix1_perX = ix1*perX
                        
                        for ixp in range(xNpMeshRes):
                            doZeroFxp = False
                            xpMeshRes = xStartMeshRes + ixp*xStepMeshRes
                            ixp0 = 0
                            if(xStepWfr > 0): ixp0 = int(trunc((xpMeshRes - xStartWfr)/xStepWfr + 1.e-09))
                            if((ixp0 < 0) or (ixp0 > xNpWfr_mi_1)):
                                doZeroFxp = True
                                #self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                #continue
                            ixp1 = ixp0 + 1
                            if(ixp1 > xNpWfr_mi_1): ixp1 = xNpWfr_mi_1
                            txp = 0
                            if(xStepWfr > 0): txp = (xpMeshRes - (xStartWfr + xStepWfr*ixp0))/xStepWfr

                            ixp0_perXp = ixp0*perXp
                            ixp1_perXp = ixp1*perXp

                            ixp0_perXp_p_ix0_perX = ixp0_perXp + ix0_perX
                            ixp1_perXp_p_ix0_perX = ixp1_perXp + ix0_perX
                            ixp0_perXp_p_ix1_perX = ixp0_perXp + ix1_perX
                            ixp1_perXp_p_ix1_perX = ixp1_perXp + ix1_perX

                            ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY = ixp0_perXp_p_ix0_perX + iyp0_perYp_p_iy0_perY
                            ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY = ixp1_perXp_p_ix0_perX + iyp0_perYp_p_iy0_perY
                            ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY = ixp0_perXp_p_ix1_perX + iyp0_perYp_p_iy0_perY
                            ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY = ixp0_perXp_p_ix0_perX + iyp1_perYp_p_iy0_perY
                            ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY = ixp0_perXp_p_ix0_perX + iyp0_perYp_p_iy1_perY

                            ixp1_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY = ixp1_perXp_p_ix1_perX + iyp0_perYp_p_iy0_perY
                            ixp1_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY = ixp1_perXp_p_ix0_perX + iyp1_perYp_p_iy0_perY
                            ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY = ixp1_perXp_p_ix0_perX + iyp0_perYp_p_iy1_perY
                            ixp0_perXp_p_ix1_perX_p_iyp1_perYp_p_iy0_perY = ixp0_perXp_p_ix1_perX + iyp1_perYp_p_iy0_perY
                            ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy1_perY = ixp0_perXp_p_ix1_perX + iyp0_perYp_p_iy1_perY
                            ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy1_perY = ixp0_perXp_p_ix0_perX + iyp1_perYp_p_iy1_perY
                            
                            for ie in range(eNpMeshRes):
                                doZeroFe = False
                                eMeshRes = eStartMeshRes + ie*eStepMeshRes
                                ie0 = 0
                                if(eStepWfr > 0): ie0 = int(trunc((eMeshRes - eStartWfr)/eStepWfr + 1.e-09))
                                if((ie0 < 0) or (ie0 > eNpWfr_mi_1)):
                                    doZeroFe = True
                                    #self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                    #continue
                                ie1 = ie0 + 1
                                if(ie1 > eNpWfr_mi_1): ie1 = eNpWfr_mi_1
                                te = 0
                                if(eStepWfr > 0): te = (eMeshRes - (eStartWfr + eStepWfr*ie0))/eStepWfr

                                ie0_perE = ie0*perE
                                ie1_perE = ie1*perE
                                
                                for iep in range(eNpMeshRes):
                                    doZeroFep = False
                                    epMeshRes = eStartMeshRes + iep*eStepMeshRes
                                    iep0 = 0
                                    if(eStepWfr > 0): iep0 = int(trunc((epMeshRes - eStartWfr)/eStepWfr + 1.e-09))
                                    if((iep0 < 0) or (iep0 > eNpWfr_mi_1)):
                                        doZeroFep = True
                                        #self.arS[ir] = self.arS[ir]*_iter/(_iter + 1); ir += 1
                                        #continue
                                    iep1 = iep0 + 1
                                    if(iep1 > eNpWfr_mi_1): iep1 = eNpWfr_mi_1
                                    tep = 0
                                    if(eStepWfr > 0): tep = (epMeshRes - (eStartWfr + eStepWfr*iep0))/eStepWfr

                                    iep0_perEp = iep0*perEp #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)
                                    iep1_perEp = iep1*perEp

                                    fInterp = 0
                                    if(not(doZeroFy or doZeroFyp or doZeroFx or doZeroFxp or doZeroFe or doZeroFep)):

                                        iep0_perEp_p_ie0_perE = iep0_perEp + ie0_perE #OC04052018 
                                        iep1_perEp_p_ie0_perE = iep1_perEp + ie0_perE
                                        iep0_perEp_p_ie1_perE = iep0_perEp + ie1_perE
                                        iep1_perEp_p_ie1_perE = iep1_perEp + ie1_perE

                                        ofst000000 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY #OC04052018 
                                        ofst100000 = iOfstSt + iep1_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst010000 = iOfstSt + iep0_perEp_p_ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst001000 = iOfstSt + iep0_perEp_p_ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst000100 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst000010 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY
                                        ofst000001 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY

                                        ofst110000 = iOfstSt + iep1_perEp_p_ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst101000 = iOfstSt + iep1_perEp_p_ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst100100 = iOfstSt + iep1_perEp_p_ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst100010 = iOfstSt + iep1_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY
                                        ofst100001 = iOfstSt + iep1_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY

                                        ofst011000 = iOfstSt + iep0_perEp_p_ie1_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst010100 = iOfstSt + iep0_perEp_p_ie1_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst010010 = iOfstSt + iep0_perEp_p_ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY
                                        ofst010001 = iOfstSt + iep0_perEp_p_ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY
                                   
                                        ofst001100 = iOfstSt + iep0_perEp_p_ie0_perE + ixp1_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY
                                        ofst001010 = iOfstSt + iep0_perEp_p_ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY
                                        ofst001001 = iOfstSt + iep0_perEp_p_ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY

                                        ofst000110 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp1_perYp_p_iy0_perY
                                        ofst000101 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy1_perY

                                        ofst000011 = iOfstSt + iep0_perEp_p_ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy1_perY

                                        #a000000 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f100000 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f010000 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f001000 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f000100 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f000010 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        #f000001 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]

                                        #f110000 = _more_stokes.arS[iOfstSt + iep1 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f101000 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f100100 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f100010 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        #f100001 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]

                                        #f011000 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f010100 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f010010 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        #f010001 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]
                                   
                                        #f001100 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        #f001010 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        #f001001 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]

                                        #f000110 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp1_perYp_p_iy0_perY]
                                        #f000101 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy1_perY]

                                        #f000011 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy1_perY]

                                        for ii in range(2): #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)

                                            a000000 = _more_stokes.arS[ofst000000]; ofst000000 += 1
                                            f100000 = _more_stokes.arS[ofst100000]; ofst100000 += 1
                                            f010000 = _more_stokes.arS[ofst010000]; ofst010000 += 1
                                            f001000 = _more_stokes.arS[ofst001000]; ofst001000 += 1
                                            f000100 = _more_stokes.arS[ofst000100]; ofst000100 += 1
                                            f000010 = _more_stokes.arS[ofst000010]; ofst000010 += 1
                                            f000001 = _more_stokes.arS[ofst000001]; ofst000001 += 1

                                            f110000 = _more_stokes.arS[ofst110000]; ofst110000 += 1
                                            f101000 = _more_stokes.arS[ofst101000]; ofst101000 += 1
                                            f100100 = _more_stokes.arS[ofst100100]; ofst100100 += 1
                                            f100010 = _more_stokes.arS[ofst100010]; ofst100010 += 1
                                            f100001 = _more_stokes.arS[ofst100001]; ofst100001 += 1

                                            f011000 = _more_stokes.arS[ofst011000]; ofst011000 += 1
                                            f010100 = _more_stokes.arS[ofst010100]; ofst010100 += 1
                                            f010010 = _more_stokes.arS[ofst010010]; ofst010010 += 1
                                            f010001 = _more_stokes.arS[ofst010001]; ofst010001 += 1

                                            f001100 = _more_stokes.arS[ofst001100]; ofst001100 += 1
                                            f001010 = _more_stokes.arS[ofst001010]; ofst001010 += 1
                                            f001001 = _more_stokes.arS[ofst001001]; ofst001001 += 1

                                            f000110 = _more_stokes.arS[ofst000110]; ofst000110 += 1
                                            f000101 = _more_stokes.arS[ofst000101]; ofst000101 += 1

                                            f000011 = _more_stokes.arS[ofst000011]; ofst000011 += 1

                                            a100000 = f100000 - a000000 
                                            a010000 = f010000 - a000000
                                            a001000 = f001000 - a000000
                                            a000100 = f000100 - a000000
                                            a000010 = f000010 - a000000
                                            a000001 = f000001 - a000000
                                            a110000 = a000000 - f010000 - f100000 + f110000
                                            a101000 = a000000 - f001000 - f100000 + f101000
                                            a100100 = a000000 - f000100 - f100000 + f100100
                                            a100010 = a000000 - f000010 - f100000 + f100010
                                            a100001 = a000000 - f000001 - f100000 + f100001
                                            a011000 = a000000 - f001000 - f010000 + f011000
                                            a010100 = a000000 - f000100 - f010000 + f010100
                                            a010010 = a000000 - f000010 - f010000 + f010010
                                            a010001 = a000000 - f000001 - f010000 + f010001
                                            a001100 = a000000 - f000100 - f001000 + f001100
                                            a001010 = a000000 - f000010 - f001000 + f001010
                                            a001001 = a000000 - f000001 - f001000 + f001001
                                            a000110 = a000000 - f000010 - f000100 + f000110
                                            a000101 = a000000 - f000001 - f000100 + f000101
                                            a000011 = a000000 - f000001 - f000010 + f000011

                                            fInterp = (a100000 + a110000*te + a101000*txp + a100100*tx + a100010*typ + a100001*ty)*tep
                                            fInterp += (a010000 + a011000*txp + a010100*tx + a010010*typ + a010001*ty)*te
                                            fInterp += (a001000 + a001100*tx + a001010*typ + a001001*ty)*txp
                                            fInterp += (a000100 + a000110*typ + a000101*ty)*tx + (a000010 + a000011*ty)*typ + a000001*ty + a000000

                                            #DEBUG
                                            #if(fInterp != 0): print(fInterp)

                                            #self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)

                                            #OC13112020
                                            if(_sum): self.arS[ir] += _mult*fInterp
                                            else: self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/iter_p_1

                                            #DEBUG
                                            #print('   ir=',ir, 'arS[ir]=', self.arS[ir])
                                           
                                            ir += 1
                                    else: #OC04052018
                                        #self.arS[ir] = self.arS[ir]*iter_d_iter_p_1; ir += 1
                                        #self.arS[ir] = self.arS[ir]*iter_d_iter_p_1
                                        self.arS[ir] *= iter_d_iter_p_1; ir += 1
                                        self.arS[ir] *= iter_d_iter_p_1; ir += 1

                                    #OC04052018 (commented-out)
                                    #self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)
                                    #ir += 1
            iOfstSt += nRadWfr        

    def to_int(self, _pol=6):
        """Calculates / "extracts" intensity at a given polarization from the Stokes components
        :param _pol: polarization component to extract: 
            0- Linear Horizontal; 
            1- Linear Vertical; 
            2- Linear 45 degrees; 
            3- Linear 135 degrees; 
            4- Circular Right; 
            5- Circular Left; 
            6- Total
        :return: 1D array with (C-aligned) resulting intensity data
        """

        resArI = None
        if(self.mutual == 0):
            nPer = self.mesh.ne*self.mesh.nx*self.mesh.ny

            lenArS = len(self.arS) #OC23072021
            if(lenArS == nPer):
                resArI = self.arS
                if(_pol != 6): print('WARNING: requested polarization component may not be estimated correctly from given Stokes parameters data')
            elif(lenArS == 4*nPer):
                resArI = array(self.numTypeStokes, [0]*nPer)
                nPer2 = 2*nPer
                nPer3 = 3*nPer
                for i in range(nPer):
                    s0 = self.arS[i]
                    s1 = self.arS[i + nPer]
                    s2 = self.arS[i + nPer2]
                    s3 = self.arS[i + nPer3]
                    resArI[i] = 0
                    if(_pol == 0): resArI[i] = 0.5*(s0 + s1) #LH
                    elif(_pol == 1): resArI[i] = 0.5*(s0 - s1) #LV
                    elif(_pol == 2): resArI[i] = 0.5*(s0 + s2) #L45
                    elif(_pol == 3): resArI[i] = 0.5*(s0 - s2) #L135
                    elif(_pol == 4): resArI[i] = 0.5*(s0 + s3) #CR (?)
                    elif(_pol == 5): resArI[i] = 0.5*(s0 - s3) #CL (?)
                    elif(_pol == 6): resArI[i] = s0 #Total
            else: raise Exception("Not all Stokes components are available in that structure")  #OC23072021
            
        #else: #to add the case of mutual intensity (what to extract: normal or mutual intensity at a given polarization?)
        return resArI

    def to_deg_coh(self, _rel_zer_tol=1.e-04, _rot=True): #05052018
        """Calculates / "extracts" Degree of Coherence from the Mutual Intensity (first Stokes component, s0)
        :param _rel_zer_tol: relative zero tolerance to use at normalizing (dividing) by the intensity
        :param _rot: rotate or not the degree of coherence data
        :return: 1D array with (C-aligned) resulting degree of coherence data
        """

        if(self.mutual <= 0): raise Exception("Calculation of Degree of Coherence can not be done from regular Intensity; Mutual Intensity is required.")

        origNe = self.mesh.ne
        origNx = self.mesh.nx
        origNy = self.mesh.ny

        origPerEp = 2
        origPerE = origPerEp*origNe
        origPerXp = origPerE*origNe
        origPerX = origPerXp*origNx
        origPerYp = origPerX*origNx
        origPerY = origPerYp*origNy

        #orig_iec = int(round(0.5*origNe)) #NOTE: int(round(0.5*1)) produces different result by Py 2.7 and Py 3.6 !!??:)
        #orig_ixc = int(round(0.5*origNx))
        #orig_iyc = int(round(0.5*origNy))
        #OC14112019
        orig_iec = 0 if(origNe <= 1) else int(round(0.5*(origNe + 1e-12)))
        orig_ixc = 0 if(origNx <= 1) else int(round(0.5*(origNx + 1e-12)))
        orig_iyc = 0 if(origNy <= 1) else int(round(0.5*(origNy + 1e-12)))
        
        ofstIc = orig_iec*origPerEp + orig_iec*origPerE + orig_ixc*origPerXp + orig_ixc*origPerX + orig_iyc*origPerYp + orig_iyc*origPerY
        absZerTolI = abs(self.arS[ofstIc])*_rel_zer_tol

        #DEBUG
        #print('   _rel_zer_tol=', _rel_zer_tol, ' absZerTolI=', absZerTolI)
        #DEBUG
        #dcWasNotZero = False
        #normIntWasZeroWhere_dcWasNotZero = False

        auxNp = origNe*origNx*origNy
        auxNp *= auxNp
        resDegCohNonRot = array('f', [0]*auxNp)

        ie2_0_origPerEp_p_ie1_0_origPerE = 0
        ie1_0_origPerEp_p_ie1_0_origPerE = 0
        ie2_0_origPerEp_p_ie2_0_origPerE = 0

        ix2_0_origPerXp_p_ix1_0_origPerX = 0
        ix1_0_origPerXp_p_ix1_0_origPerX = 0
        ix2_0_origPerXp_p_ix2_0_origPerX = 0

        iy2_0_origPerYp_p_iy1_0_origPerY = 0
        iy1_0_origPerYp_p_iy1_0_origPerY = 0
        iy2_0_origPerYp_p_iy2_0_origPerY = 0

        resOfst = 0
        for iy in range(origNy):
            for iyp in range(origNy):
                if(origNy > 1):
                    iy1_0_origPerY = iy*origPerY
                    iy2_0_origPerYp = iyp*origPerYp
                    iy2_0_origPerYp_p_iy1_0_origPerY = iy2_0_origPerYp + iy1_0_origPerY
                    iy1_0_origPerYp_p_iy1_0_origPerY = iy*origPerYp + iy1_0_origPerY
                    iy2_0_origPerYp_p_iy2_0_origPerY = iy2_0_origPerYp + iyp*origPerY

                for ix in range(origNx):
                    for ixp in range(origNx):
                        if(origNx > 1):
                            ix1_0_origPerX = ix*origPerX
                            ix2_0_origPerXp = ixp*origPerXp
                            ix2_0_origPerXp_p_ix1_0_origPerX = ix2_0_origPerXp + ix1_0_origPerX
                            ix1_0_origPerXp_p_ix1_0_origPerX = ix*origPerXp + ix1_0_origPerX
                            ix2_0_origPerXp_p_ix2_0_origPerX = ix2_0_origPerXp + ixp*origPerX

                        if(origNe > 1):
                            for ie in range(origNe):
                                for iep in range(origNe):
                                    ie1_0_origPerE = ie*origPerE
                                    ie2_0_origPerEp = iep*origPerEp
                                    ie2_0_origPerEp_p_ie1_0_origPerE = ie2_0_origPerEp + ie1_0_origPerE
                                    ie1_0_origPerEp_p_ie1_0_origPerE = ie*origPerEp + ie1_0_origPerE
                                    ie2_0_origPerEp_p_ie2_0_origPerE = ie2_0_origPerEp + iep*origPerE

                                    ofstMI = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                    ofstI1 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                    ofstI2 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY

                                    reMI = self.arS[ofstMI]; imMI = self.arS[ofstMI + 1]
                                    absMI = sqrt(reMI*reMI + imMI*imMI)
                                    reI1 = self.arS[ofstI1]#; imI1 = self.arS[ofstI1 + 1]
                                    reI2 = self.arS[ofstI2]#; imI2 = self.arS[ofstI2 + 1]

                                    denom = sqrt(reI1*reI2) + absZerTolI #OC31072018
                                    if(denom == 0): resDegCohNonRot[resOfst] = 0
                                    else: resDegCohNonRot[resOfst] = absMI/denom

                                    #resDegCohNonRot[resOfst] = absMI/(sqrt(reI1*reI2) + absZerTolI)
                                    
                                    resOfst += 1
                        else:
                            ofstMI = ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                            
                            ofstI1 = ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                            ofstI2 = ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY

                            reMI = float(self.arS[ofstMI]); imMI = float(self.arS[ofstMI + 1]) #OC01122021
                            #reMI = self.arS[ofstMI]; imMI = self.arS[ofstMI + 1]

                            absMI = sqrt(reMI*reMI + imMI*imMI)
                            reI1 = float(self.arS[ofstI1]) #OC01122021
                            reI2 = float(self.arS[ofstI2]) #OC01122021
                            #reI1 = self.arS[ofstI1]#; imI1 = self.arS[ofstI1 + 1]
                            #reI2 = self.arS[ofstI2]#; imI2 = self.arS[ofstI2 + 1]

                            denom = sqrt(reI1*reI2) + absZerTolI #OC31072018
                            if(denom == 0): resDegCohNonRot[resOfst] = 0
                            else: resDegCohNonRot[resOfst] = absMI/denom

                            #DEBUG
                            #if(absMI > 0):
                            #    dcWasNotZero = True
                            #    if(denom == 0): normIntWasZeroWhere_dcWasNotZero = True
                          
                            #resDegCohNonRot[resOfst] = absMI/(sqrt(reI1*reI2) + absZerTolI)
                            
                            resOfst += 1

        #DEBUG
        #if(dcWasNotZero):
        #    print('DC was NOT ZERO at some points')
        #    if(normIntWasZeroWhere_dcWasNotZero): print('...BUT Normal intensity WAS ZERO there')
        #else: print('DC was ZERO at all points')
        #DEBUG
        #rotDC_WasNotZero = False
                            
        if(not _rot): return resDegCohNonRot
        #DEBUG
        #return resDegCohNonRot

        origEstart = self.mesh.eStart
        origEfin = self.mesh.eFin
        origXstart = self.mesh.xStart
        origXfin = self.mesh.xFin
        origYstart = self.mesh.yStart
        origYfin = self.mesh.yFin

        origEstep = 0; origEstepInv = 0
        if(origNe > 1): 
            origEstep = (origEfin - origEstart)/(origNe - 1)
            origEstepInv = 1/origEstep
        origXstep = 0; origXstepInv = 0
        if(origNx > 1): 
            origXstep = (origXfin - origXstart)/(origNx - 1)
            origXstepInv = 1/origXstep
        origYstep = 0; origYstepInv = 0
        if(origNy > 1): 
            origYstep = (origYfin - origYstart)/(origNy - 1)
            origYstepInv = 1/origYstep

        origNe_m_1 = origNe - 1; origNe_m_2 = origNe - 2
        origNx_m_1 = origNx - 1; origNx_m_2 = origNx - 2
        origNy_m_1 = origNy - 1; origNy_m_2 = origNy - 2

        resNe = origNe
        resNx = origNx
        resNy = origNy
        
        #resEstart = origEstart
        #resEfin = origEfin
        #resXstart = origXstart
        #resXfin = origXfin
        #resYstart = origYstart
        #resYfin = origYfin
        
        #OC12112019
        resE_start = origEstart
        resE_fin = origEfin
        resEp_start = 0.5*(origEstart - origEfin)
        #resEp_fin = 0.5*(origEfin - origEstart)
        resX_start = origXstart
        resX_fin = origXfin
        resXp_start = 0.5*(origXstart - origXfin)
        #resXp_fin = 0.5*(origXfin - origXstart)
        resY_start = origYstart
        resY_fin = origYfin
        resYp_start = 0.5*(origYstart - origYfin)
        #resYp_fin = 0.5*(origYfin - origYstart)

        #sqrt2 = sqrt(2)
        #if(resNe > 1): 
        #    resNe = 2*resNe - 1
        #    ec = 0.5*(resEstart + resEfin)
        #    resEstart = ec - sqrt2*(ec - resEstart)
        #    resEfin = ec + sqrt2*(resEfin - ec)
        #if(resNx > 1): 
        #    resNx = 2*resNx - 1
        #    xc = 0.5*(resXstart + resXfin)
        #    resXstart = xc - sqrt2*(xc - resXstart)
        #    resXfin = xc + sqrt2*(resXfin - xc)
        #if(resNy > 1): 
        #    resNy = 2*resNy - 1
        #    yc = 0.5*(resYstart + resYfin)
        #    resYstart = yc - sqrt2*(yc - resYstart)
        #    resYfin = yc + sqrt2*(resYfin - yc)

        resNp = resNe*resNx*resNy
        resNp *= resNp
        resDegCoh = array('f', [0]*resNp)

        #resEstep = 0
        #if(resNe > 1): resEstep = (resEfin - resEstart)/(resNe - 1)
        #resXstep = 0
        #if(resNx > 1): resXstep = (resXfin - resXstart)/(resNx - 1)
        #resYstep = 0
        #if(resNy > 1): resYstep = (resYfin - resYstart)/(resNy - 1)

        #OC12112019
        resEstep = 0
        if(resNe > 1): resEstep = (resE_fin - resE_start)/(resNe - 1)
        resXstep = 0
        if(resNx > 1): resXstep = (resX_fin - resX_start)/(resNx - 1)
        resYstep = 0
        if(resNy > 1): resYstep = (resY_fin - resY_start)/(resNy - 1)

        perE = origNe
        perXp = perE*origNe
        perX = perXp*origNx
        perYp = perX*origNx
        perY = perYp*origNy

        resX1 = 0; resX2 = 0 #just declaring these vars
        resE1 = 0; resE2 = 0 

        ie2_0_p_ie1_0_perE = 0
        ie2_1_p_ie1_0_perE = 0
        ie2_0_p_ie1_1_perE = 0
        ie2_1_p_ie1_1_perE = 0

        ix2_0_perXp_p_ix1_0_perX = 0
        ix2_1_perXp_p_ix1_0_perX = 0
        ix2_0_perXp_p_ix1_1_perX = 0
        ix2_1_perXp_p_ix1_1_perX = 0

        iy2_0_perYp_p_iy1_0_perY = 0
        iy2_1_perYp_p_iy1_0_perY = 0
        iy2_0_perYp_p_iy1_1_perY = 0
        iy2_1_perYp_p_iy1_1_perY = 0

        resOfst = 0
        #resYp = resYstart #resY = resYstart
        resYp = resYp_start #OC12112019
        for iyp in range(resNy): #for iy in range(resNy):
            #resY = resYstart #resYp = resYstart
            resY = resY_start #OC12112019
            for iy in range(resNy): #for iyp in range(resNy):
                isInRangeY = True
                if(resNy > 1):
                    resY1 = resY + resYp
                    resY2 = resY - resYp
                    #resY1 = resY - resYp
                    #resY2 = resY + resYp
                    if((resY1 < origYstart) or (resY1 > origYfin) or (resY2 < origYstart) or (resY2 > origYfin)): isInRangeY = False

                iy1_0 = 0; iy1_1 = 0; ry1 = 0
                iy2_0 = 0; iy2_1 = 0; ry2 = 0
                if((resNy > 1) and isInRangeY):
                    iy1_0 = int(trunc((resY1 - origYstart)*origYstepInv))
                    if(iy1_0 >= origNy_m_1): iy1_0 = origNy_m_2
                    ry1 = (resY1 - (origYstart + iy1_0*origYstep))*origYstepInv

                    iy2_0 = int(trunc((resY2 - origYstart)*origYstepInv))
                    if(iy2_0 >= origNy_m_1): iy2_0 = origNy_m_2
                    ry2 = (resY2 - (origYstart + iy2_0*origYstep))*origYstepInv

                    iy1_1 = iy1_0 + 1
                    if(iy1_1 > origNy_m_1): iy1_1 = origNy_m_1
                    iy2_1 = iy2_0 + 1
                    if(iy2_1 > origNy_m_1): iy2_1 = origNy_m_1

                    iy2_0_perYp = iy2_0*perYp
                    iy2_1_perYp = iy2_1*perYp
                    iy1_0_perY = iy1_0*perY
                    iy1_1_perY = iy1_1*perY
                    iy2_0_perYp_p_iy1_0_perY = iy2_0_perYp + iy1_0_perY
                    iy2_1_perYp_p_iy1_0_perY = iy2_1_perYp + iy1_0_perY
                    iy2_0_perYp_p_iy1_1_perY = iy2_0_perYp + iy1_1_perY
                    iy2_1_perYp_p_iy1_1_perY = iy2_1_perYp + iy1_1_perY

                #resXp = resXstart #resX = resXstart
                resXp = resXp_start #OC12112019
                for ixp in range(resNx): #for ix in range(resNx):
                    #resX = resXstart #resXp = resXstart
                    resX = resX_start #OC12112019
                    for ix in range(resNx): #for ixp in range(resNx):
                        isInRangeX = True
                        if(isInRangeY):
                            if(resNx > 1):
                                resX1 = resX + resXp
                                resX2 = resX - resXp
                                #resX1 = resX - resXp
                                #resX2 = resX + resXp
                                if((resX1 < origXstart) or (resX1 > origXfin) or (resX2 < origXstart) or (resX2 > origXfin)): isInRangeX = False
                        else:
                            isInRangeX = False

                        ix1_0 = 0; ix1_1 = 0; rx1 = 0
                        ix2_0 = 0; ix2_1 = 0; rx2 = 0
                        if((resNx > 1) and isInRangeX):
                            ix1_0 = int(trunc((resX1 - origXstart)*origXstepInv))
                            if(ix1_0 >= origNx_m_1): ix1_0 = origNx_m_2
                            rx1 = (resX1 - (origXstart + ix1_0*origXstep))*origXstepInv

                            ix2_0 = int(trunc((resX2 - origXstart)*origXstepInv))
                            if(ix2_0 >= origNx_m_1): ix2_0 = origNx_m_2
                            rx2 = (resX2 - (origXstart + ix2_0*origXstep))*origXstepInv

                            ix1_1 = ix1_0 + 1
                            if(ix1_1 > origNx_m_1): ix1_1 = origNx_m_1
                            ix2_1 = ix2_0 + 1
                            if(ix2_1 > origNx_m_1): ix2_1 = origNx_m_1

                            ix2_0_perXp = ix2_0*perXp
                            ix2_1_perXp = ix2_1*perXp
                            ix1_0_perX = ix1_0*perX
                            ix1_1_perX = ix1_1*perX
                            ix2_0_perXp_p_ix1_0_perX = ix2_0_perXp + ix1_0_perX
                            ix2_1_perXp_p_ix1_0_perX = ix2_1_perXp + ix1_0_perX
                            ix2_0_perXp_p_ix1_1_perX = ix2_0_perXp + ix1_1_perX
                            ix2_1_perXp_p_ix1_1_perX = ix2_1_perXp + ix1_1_perX

                        #resEp = resEstart #resE = resEstart
                        resEp = resEp_start #OC12112019
                        for iep in range(resNe): #for ie in range(resNe):
                            #resE = resEstart #resEp = resEstart
                            resE = resE_start #OC12112019                            
                            for ie in range(resNe): #for iep in range(resNe):
                                isInRangeE = True
                                if(isInRangeX):
                                    if(resNe > 1):
                                        resE1 = resE + resEp
                                        resE2 = resE - resEp
                                        #resE1 = resE - resEp
                                        #resE2 = resE + resEp
                                        if((resE1 < origEstart) or (resE1 > origEfin) or (resE2 < origEstart) or (resE2 > origEfin)): isInRangeE = False
                                else:
                                    isInRangeE = False

                                ie1_0 = 0; ie1_1 = 0; re1 = 0
                                ie2_0 = 0; ie2_1 = 0; re2 = 0
                                if((resNe > 1) and isInRangeE):
                                    ie1_0 = int(trunc((resE1 - origEstart)*origEstepInv))
                                    if(ie1_0 >= origNe_m_1): ie1_0 = origNe_m_2
                                    re1 = (resE1 - (origEstart + ie1_0*origEstep))*origEstepInv

                                    ie2_0 = int(trunc((resE2 - origEstart)*origEstepInv))
                                    if(ie2_0 >= origNe_m_1): ie2_0 = origNe_m_2
                                    re2 = (resE2 - (origEstart + ie2_0*origEstep))*origEstepInv

                                    ie1_1 = ie1_0 + 1
                                    if(ie1_1 > origNe_m_1): ie1_1 = origNe_m_1
                                    ie2_1 = ie2_0 + 1
                                    if(ie2_1 > origNe_m_1): ie2_1 = origNe_m_1

                                    ie1_0_perE = ie1_0*perE
                                    ie1_1_perE = ie1_1*perE
                                    ie2_0_p_ie1_0_perE = ie2_0 + ie1_0_perE
                                    ie2_1_p_ie1_0_perE = ie2_1 + ie1_0_perE
                                    ie2_0_p_ie1_1_perE = ie2_0 + ie1_1_perE
                                    ie2_1_p_ie1_1_perE = ie2_1 + ie1_1_perE

                                if(isInRangeE):
                                    a000000 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f100000 = resDegCohNonRot[ie2_1_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f010000 = resDegCohNonRot[ie2_0_p_ie1_1_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f001000 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_1_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f000100 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_1_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f000010 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_1_perYp_p_iy1_0_perY]
                                    f000001 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_1_perY]

                                    f110000 = resDegCohNonRot[ie2_1_p_ie1_1_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f101000 = resDegCohNonRot[ie2_1_p_ie1_0_perE + ix2_1_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f100100 = resDegCohNonRot[ie2_1_p_ie1_0_perE + ix2_0_perXp_p_ix1_1_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f100010 = resDegCohNonRot[ie2_1_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_1_perYp_p_iy1_0_perY]
                                    f100001 = resDegCohNonRot[ie2_1_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_1_perY]

                                    f011000 = resDegCohNonRot[ie2_0_p_ie1_1_perE + ix2_1_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f010100 = resDegCohNonRot[ie2_0_p_ie1_1_perE + ix2_0_perXp_p_ix1_1_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f010010 = resDegCohNonRot[ie2_0_p_ie1_1_perE + ix2_0_perXp_p_ix1_0_perX + iy2_1_perYp_p_iy1_0_perY]
                                    f010001 = resDegCohNonRot[ie2_0_p_ie1_1_perE + ix2_0_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_1_perY]

                                    f001100 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_1_perXp_p_ix1_1_perX + iy2_0_perYp_p_iy1_0_perY]
                                    f001010 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_1_perXp_p_ix1_0_perX + iy2_1_perYp_p_iy1_0_perY]
                                    f001001 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_1_perXp_p_ix1_0_perX + iy2_0_perYp_p_iy1_1_perY]

                                    f000110 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_1_perX + iy2_1_perYp_p_iy1_0_perY]
                                    f000101 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_1_perX + iy2_0_perYp_p_iy1_1_perY]

                                    f000011 = resDegCohNonRot[ie2_0_p_ie1_0_perE + ix2_0_perXp_p_ix1_0_perX + iy2_1_perYp_p_iy1_1_perY]

                                    a100000 = f100000 - a000000
                                    a010000 = f010000 - a000000
                                    a001000 = f001000 - a000000
                                    a000100 = f000100 - a000000
                                    a000010 = f000010 - a000000
                                    a000001 = f000001 - a000000
                                    a110000 = a000000 - f010000 - f100000 + f110000
                                    a101000 = a000000 - f001000 - f100000 + f101000
                                    a100100 = a000000 - f000100 - f100000 + f100100
                                    a100010 = a000000 - f000010 - f100000 + f100010
                                    a100001 = a000000 - f000001 - f100000 + f100001
                                    a011000 = a000000 - f001000 - f010000 + f011000
                                    a010100 = a000000 - f000100 - f010000 + f010100
                                    a010010 = a000000 - f000010 - f010000 + f010010
                                    a010001 = a000000 - f000001 - f010000 + f010001
                                    a001100 = a000000 - f000100 - f001000 + f001100
                                    a001010 = a000000 - f000010 - f001000 + f001010
                                    a001001 = a000000 - f000001 - f001000 + f001001
                                    a000110 = a000000 - f000010 - f000100 + f000110
                                    a000101 = a000000 - f000001 - f000100 + f000101
                                    a000011 = a000000 - f000001 - f000010 + f000011

                                    rotDegCoh = (a100000 + a110000*re1 + a101000*rx2 + a100100*rx1 + a100010*ry2 + a100001*ry1)*re2
                                    rotDegCoh += (a010000 + a011000*rx2 + a010100*rx1 + a010010*ry2 + a010001*ry1)*re1
                                    rotDegCoh += (a001000 + a001100*rx1 + a001010*ry2 + a001001*ry1)*rx2
                                    rotDegCoh += (a000100 + a000110*ry2 + a000101*ry1)*rx1 + (a000010 + a000011*ry1)*ry2 + a000001*ry1 + a000000
                                    #OC05052018: Note that this does not include all terms of multi-dim "bi-linear" interpolation!?

                                    resDegCoh[resOfst] = rotDegCoh

                                    #DEBUG
                                    #if(rotDegCoh > 0): rotDC_WasNotZero = True
                                    
                                else:
                                    resDegCoh[resOfst] = 0

                                resOfst += 1
                                resE += resEstep #resEp += resEstep
                            resEp += resEstep #resE += resEstep
                        resX += resXstep #resXp += resXstep
                    resXp += resXstep #resX += resXstep
                resY += resYstep #resYp += resYstep
            resYp += resYstep #resY += resYstep

        #DEBUG
        #if(rotDC_WasNotZero): print('Rotated DC WAS NOT ZERO at some points')
        #else: print('Rotated DC WAS ZERO at all points')
        
        return resDegCoh

    def extr_save_pol_comp(self, _type=2, _pol=6, _fname=None, _fname_stk=None): #OC27122023
    #def extr_save_spec_comp(self, stk, _type=2, _pol=6, _fname=None, _fname_stk=None): #OC10102023
        """Extracts Flux / Intensity at different polarizations from Stokes data and saves it to a file"""
        
        arI = None
        if(_pol < 7):
            arI = self.to_int(_pol) #OC27122023
            #arI = stk.to_int(_pol)
        else:
            arI = []
            for i in range(6):
                arI.append(self.to_int(i)) #OC27122023
                #arI.append(stk.to_int(i))
            
        if(srwl_uti_proc_is_master()): #Consider removing it from here?
            sValName = 'Flux'
            sValUnitName = 'ph/s/.1%bw'
            if(_type == 2):
                sValName = 'Intensity'
                sValUnitName = 'ph/s/.1%bw/mm^2'

            if(_pol != 6): 
                if(_fname_stk is not None):
                    if(len(_fname_stk) > 0):
                        srwl_uti_save_intens_ascii(self.arS, self.mesh, _fname_stk, _n_stokes=4, _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName]) #OC27122023
                        #srwl_uti_save_intens_ascii(stk.arS, stk.mesh, _fname_stk, _n_stokes=4, _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
                
            if(_fname is not None):
                if(len(_fname) > 0):
                    #sValName = 'Flux'
                    #sValUnitName = 'ph/s/.1%bw'
                    #if(_type == 2):
                    #    sValName = 'Intensity'
                    #    sValUnitName = 'ph/s/.1%bw/mm^2'
                        
                    if(_pol < 7): #OC10102023
                        srwl_uti_save_intens_ascii(arI, self.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName]) #OC27122023
                        #srwl_uti_save_intens_ascii(arI, stk.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
                    else: #OC10102023 (adding "_0" .. "_5" before ".dat" here)
                        len_fname = len(_fname)
                        indLastDot = _fname.rfind('.')
                        nameCore = ''; sExt = ''
                        if((indLastDot >= 0) and (indLastDot < len_fname)):
                            nameCore = _fname[:indLastDot]
                            sExt = _fname[indLastDot:len_fname] #extension with '.'
                        for i in range(6):
                            fnPol = nameCore + '_' + repr(i) + sExt
                            srwl_uti_save_intens_ascii(arI[i], self.mesh, fnPol, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName]) #OC27122023
                            #srwl_uti_save_intens_ascii(arI[i], stk.mesh, fnPol, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', sValName], _arUnits=['eV', 'm', 'm', sValUnitName])
        return arI

#****************************************************************************
class SRWLWfr(object):
    """Radiation Wavefront (Electric Field)"""
    #arEx = 0 #array('f', [0]*2) #horizontal complex electric field component array; NOTE: only 'f' (float) is supported for the moment (Jan. 2011)
    #arEy = 0 #array('f', [0]*2) #vertical complex electric field component array
    #mesh = SRWLRadMesh()
    #Rx = 0 #instant wavefront radii
    #Ry = 0 
    #dRx = 0 #error of wavefront radii
    #dRy = 0
    #xc = 0 #instant transverse coordinates of wavefront instant "source center"
    #yc = 0
    #avgPhotEn = 0 #average photon energy for time-domain simulations    
    #presCA = 0 #presentation/domain: 0- coordinates, 1- angles
    #presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
    #numTypeElFld = 'f' #electric field numerical type: 'f' (float) or 'd' (double)
    #unitElFld = 1 #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)
    #partBeam = SRWLPartBeam() #particle beam source; strictly speaking, it should be just SRWLParticle; however, "multi-electron" information can appear useful for those cases when "multi-electron intensity" can be deduced from the "single-electron" one by convolution
    #arElecPropMatr = array('d', [0]*20) #effective 1st order "propagation matrix" for electron beam parameters
    #arMomX = array('d', [0]*11) #statistical moments (of Wigner distribution); to check the exact number of moments required
    #arMomY = array('d', [0]*11)
    #arWfrAuxData = array('d', [0]*30) #array of auxiliary wavefront data

    def __init__(self, _arEx=None, _arEy=None, _typeE='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0, _zStart=0, _partBeam=None):
        """
        :param _arEx: horizontal complex electric field component array; NOTE: only 'f' (float) is supported for the moment (Jan. 2011)
        :param _arEy: vertical complex electric field component array
        :param _typeE: electric field numerical type: 'f' (float) or 'd' (double)
        :param _eStart: initial value of photon energy (/time)
        :param _eFin: final value of photon energy (/time)
        :param _ne: numbers of points vs photon energy
        :param _xStart: initial value of horizontal positions
        :param _xFin: final value of horizontal positions
        :param _nx: numbers of points vs horizontal positions
        :param _yStart: initial vertical positions
        :param _yFin: final value of vertical positions
        :param _ny: numbers of points vs vertical positions
        :param _zStart: longitudinal position
        :param _partBeam: particle beam source; strictly speaking, it should be just SRWLParticle; however, "multi-electron" information can appear useful for those cases when "multi-electron intensity" can be deduced from the "single-electron" one by convolution

        Some additional parameters, that are not included in constructor arguments:
        Rx, Ry: instant wavefront radii
        dRx, dRy: error of wavefront radii
        xc, yc: transverse coordinates of wavefront instant "source center"
        avgPhotEn: average photon energy for time-domain simulations
        presCA: presentation/domain: 0- coordinates, 1- angles
        presFT: presentation/domain: 0- frequency (photon energy), 1- time
        unitElFld: electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2)
        arElecPropMatr: effective 1st order "propagation matrix" for electron beam parameters
        arMomX, arMomY: statistical moments (of Wigner distribution); to check the exact number of moments required
        arWfrAuxData: array of auxiliary wavefront data
        """
        self.arEx = _arEx
        self.arEy = _arEy
        #self.mesh = SRWLRadMesh(_eStart, _eFin, _ne, _xStart, _xFin, _nx, _yStart, _yFin, _ny)
        self.mesh = SRWLRadMesh(_eStart, _eFin, _ne, _xStart, _xFin, _nx, _yStart, _yFin, _ny, _zStart)
        self.numTypeElFld = _typeE
        self.partBeam = SRWLPartBeam() if _partBeam is None else _partBeam

        self.Rx = 0 #instant wavefront radii
        self.Ry = 0
        self.dRx = 0 #error of wavefront radii
        self.dRy = 0
        self.xc = 0 #instant transverse coordinates of wavefront instant "source center"
        self.yc = 0
        self.avgPhotEn = 0 #average photon energy for time-domain simulations
        self.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
        self.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
        self.unitElFld = 1 #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time) ?
        self.unitElFldAng = 0 #electric field units in angular representation: 0- sqrt(Wavelength[m]*Phot/s/0.1%bw/mrad^2) vs rad/Wavelength[m], 1- sqrt(Phot/s/0.1%bw/mrad^2) vs rad; [Phot/s/0.1%bw] can be replaced by [J/eV] or [W], depending on self.unitElFld, self.presFT and self.presCA
        self.arElecPropMatr = array('d', [0] * 20) #effective 1st order "propagation matrix" for electron beam parameters
        self.arMomX = array('d', [0] * 11 * _ne) #statistical moments (of Wigner distribution); to check the exact number of moments required
        self.arMomY = array('d', [0] * 11 * _ne)
        self.arWfrAuxData = array('d', [0] * 30) #array of auxiliary wavefront data

        nProd = _ne * _nx * _ny #array length to store one component of complex electric field
        EXNeeded = 0
        EYNeeded = 0
        if(_arEx == 1) and (nProd > 0):
            EXNeeded = 1
        if(_arEy == 1) and (nProd > 0):
            EYNeeded = 1
        if(EXNeeded > 0) or (EYNeeded > 0):
            self.allocate(_ne, _nx, _ny, EXNeeded, EYNeeded)

    #def allocate(self, _ne, _nx, _ny, EXNeeded=1, EYNeeded=1, typeE='f'):
    def allocate(self, _ne, _nx, _ny, _EXNeeded=1, _EYNeeded=1, _typeE='f', _backupNeeded=0): #OC141115
        """Allocate Electric Field data
        :param _ne: number of points vs photon energy / time
        :param _nx: number of points vs horizontal position / angle
        :param _ny: number of points vs vertical position / angle
        :param _EXNeeded: switch specifying whether Ex data is necessary or not (1 or 0)
        :param _EYNeeded: switch specifying whether Ey data is necessary or not (1 or 0)
        :param _typeE: numerical type of Electric Field data: float (single precision) or double ('f' or 'd'); double is not yet supported
        :param _backupNeeded: switch specifying whether backup of Electric Field data (arExAux, arEyAux) should be created or not (1 or 0)
        """
        #print('') #debugging
        #print('          (re-)allocating: old point numbers: ne=',self.mesh.ne,' nx=',self.mesh.nx,' ny=',self.mesh.ny) #,' type:',self.numTypeElFld)
        #print('                           new point numbers: ne=',_ne,' nx=',_nx,' ny=',_ny) #,' type:',typeE)
        #print('                           backupNeeded',_backupNeeded)
        
        nTot = 2*_ne*_nx*_ny #array length to store one component of complex electric field
        nMom = 11*_ne
        if _EXNeeded:
            #print('          trying to (re-)allocate Ex ... ', end='')
            #del self.arEx
            if _backupNeeded: #OC141115
                self.arExAux = self.arEx
                #print('          self.arExAux assigned') #debugging
            #else:
            #    del self.arEx
            #self.arEx = array(typeE, [0]*nTot)
            self.arEx = srwl_uti_array_alloc(_typeE, nTot)
            #print('          done')           
            if len(self.arMomX) != nMom:
                del self.arMomX
                self.arMomX = array('d', [0]*nMom)
        if _EYNeeded:
            #print('          trying to (re-)allocate Ey ... ', end='')
            #del self.arEy
            if _backupNeeded: #OC141115
                self.arEyAux = self.arEy
                #print('          self.arEyAux assigned') #debugging
            #else:
            #    del self.arEy
            #self.arEy = array(typeE, [0]*nTot)
            self.arEy = srwl_uti_array_alloc(_typeE, nTot)
            #print('          done')
            if len(self.arMomY) != nMom:
                del self.arMomY
                self.arMomY = array('d', [0]*nMom)
        self.numTypeElFld = _typeE
        self.mesh.ne = _ne
        self.mesh.nx = _nx
        self.mesh.ny = _ny
        #gc.collect() #OC14122023
        #print('                           exiting allocate')

    def delE(self, _type=0, _treatEX=1, _treatEY=1): #OC151115
        """Delete Electric Field data
        :param _type: type of data to be deleted: 0- arEx, arEy, arExAux, arEyAux; 1- arEx, arEy only; 2- arExAux, arEyAux only
        :param _treatEX: switch specifying whether Ex data should be deleted or not (1 or 0)
        :param _treatEY: switch specifying whether Ey data should be deleted or not (1 or 0)
        """
        if _treatEX:
            if((_type == 0) or (_type == 1)):
                if(self.arEx is not None):
                    del self.arEx
                    self.arEx = None
            if((_type == 0) or (_type == 2)):
                if(hasattr(self, 'arExAux')):
                    if(self.arExAux is not None):
                        del self.arExAux
                        #gc.collect() #OC13122023
                        self.arExAux = None
                        #print('          self.arExAux deleted') #debugging
        if _treatEY:
            if((_type == 0) or (_type == 1)):
                if(self.arEy is not None):
                    del self.arEy
                    self.arEy = None
            if((_type == 0) or (_type == 2)):
                if(hasattr(self, 'arEyAux')):
                    if(self.arEyAux is not None):
                        del self.arEyAux
                        #gc.collect() #OC13122023
                        self.arEyAux = None
                        #print('          self.arEyAux deleted') #debugging


    def addE(self, _wfr, _meth=0):
        """Add Another Electric Field Wavefront
        :param _wfr: wavefront to be added
        :param _meth: method of adding the wavefront _wfr:
        0- simple addition assuming _wfr to have same mesh as this wavefront
        1- add using bilinear interpolation (taking into account meshes of the two wavefronts)
        2- add using bi-quadratic interpolation (taking into account meshes of the two wavefronts)
        3- add using bi-cubic interpolation (taking into account meshes of the two wavefronts)
        """
        if(_meth == 0):
            if((self.mesh.ne != _wfr.mesh.ne) or (self.mesh.nx != _wfr.mesh.nx) or (self.mesh.ny != _wfr.mesh.ny)):
                 raise Exception("Electric Field addition can not be performed by this method because of unequal sizes of the two Wavefronts") 
            nTot = 2*self.mesh.ne*self.mesh.nx*self.mesh.ny
            #test:
            #aux = 0
            wfr_arEx = _wfr.arEx
            wfr_arEy = _wfr.arEy
            
            for i in range(nTot):
                #for some reason, this increases memory requirements in Py:
                self.arEx[i] += wfr_arEx[i] 
                self.arEy[i] += wfr_arEy[i]
                
        elif(_meth == 1):
            #to implement
            raise Exception("This Electric Field addition method is not implemented yet")
        
        elif(_meth == 2):
            #to implement
            raise Exception("This Electric Field addition method is not implemented yet")
            
        elif(_meth == 3):
            #to implement
            raise Exception("This Electric Field addition method is not implemented yet")
    
    def copy_comp(self, _stokes): 
        """Copy compenents of Electric Field to Stokes structure"""
        if(isinstance(_stokes, SRWLStokes) == False):
            raise Exception("Incorrect Stokes parameters object submitted") 
        nTot = self.mesh.ne*self.mesh.nx*self.mesh.ny
        nTotSt = nTot*4
        nTot2 = nTot*2
        nTot3 = nTot*3
        if(_stokes.arS is not None):
            if(len(_stokes.arS) < nTotSt):
                _stokes.arS = array('f', [0]*nTotSt)
        else:
            _stokes.arS = array('f', [0]*nTotSt)           
        for i in range(nTot):
            i2 = i*2
            i2p1 = i2 + 1
            #reEx = self.arEx[i2]
            #imEx = self.arEx[i2p1]
            #reEy = self.arEy[i2]
            #imEy = self.arEy[i2p1]
            _stokes.arS[i] = self.arEx[i2] #reEx
            _stokes.arS[i + nTot] = self.arEx[i2p1] #imEx
            _stokes.arS[i + nTot2] = self.arEy[i2] #reEy
            _stokes.arS[i + nTot3] = self.arEy[i2p1] #imEy
        _stokes.mesh.set_from_other(self.mesh)

    #def resize_mesh(self, _mesh): #Implemented in C++, see srwl.ResizeElecFieldMesh()
    #    """Resizes the Electric Field according to given input mesh params (_mesh)"""
    #    if(_mesh.is_equal(self.mesh)): return
    #    if((self.mesh.nx > 1) and (self.mesh.ny > 1) and (_mesh.nx > 1) and (_mesh.ny > 1)):
    #        xRange = self.mesh.xFin - self.mesh.xStart
    #        xRangeFin = _mesh.xFin - _mesh.xStart
    #        xRangeFact = xRangeFin/xRange
    #        xStep = xRange/(self.mesh.nx - 1)
    #        xStepFin = xRangeFin/(_mesh.nx - 1)
    #        xResolFact = xStep/xStepFin
    #        xcFin = 0.5*(_mesh.xStart + _mesh.xFin)
    #        xCenFact = (xcFin - self.mesh.xStart)/xRange
    #        yRange = self.mesh.yFin - self.mesh.yStart
    #        yRangeFin = _mesh.yFin - _mesh.yStart
    #        yRangeFact = yRangeFin/yRange
    #        yStep = yRange/(self.mesh.ny - 1)
    #        yStepFin = yRangeFin/(_mesh.ny - 1)
    #        yResolFact = yStep/yStepFin
    #        ycFin = 0.5*(_mesh.yStart + _mesh.yFin)
    #        yCenFact = (ycFin - self.mesh.yStart)/yRange       
    #        #DEBUG
    #        #print('')
    #        #print(xStepFin, yStepFin)
    #        #print(xcFin, self.mesh.xStart, xRange)
    #        #print(ycFin, self.mesh.yStart, yRange)
    #        #print('Resizing Params:', xRangeFact, xResolFact, yRangeFact, yResolFact, xCenFact, yCenFact)
    #        #END DEBUG
    #        srwl.ResizeElecField(self, 'c', [0, xRangeFact, xResolFact, yRangeFact, yResolFact, xCenFact, yCenFact])
    
    def calc_stokes(self, _stokes, _n_stokes_comp=4, _rx_avg=0, _ry_avg=0, _xc_avg=0, _yc_avg=0): #OC21052020
    #def calc_stokes(self, _stokes, _n_stokes_comp=4): #OC04052018
    #def calc_stokes(self, _stokes):
        """Calculate Stokes parameters from Electric Field"""
        if(_stokes.mutual <= 0):
            nTot = self.mesh.ne*self.mesh.nx*self.mesh.ny
            #if(type(_stokes).__name__ != 'SRWLStokes')):
            if(isinstance(_stokes, SRWLStokes) == False):
                raise Exception("Incorrect Stokes parameters object submitted") 
            nTotSt = nTot*_n_stokes_comp #OC18072021
            #nTotSt = nTot*4

            nTot2 = nTot*2
            nTot3 = nTot*3
            if(_stokes.arS is not None):
                if(len(_stokes.arS) < nTotSt):
                    _stokes.arS = array('f', [0]*nTotSt)
            else:
                _stokes.arS = array('f', [0]*nTotSt)           
            for i in range(nTot):
                i2 = i*2
                i2p1 = i2 + 1
                reEx = self.arEx[i2]
                imEx = self.arEx[i2p1]
                reEy = self.arEy[i2]
                imEy = self.arEy[i2p1]
                intLinX = reEx*reEx + imEx*imEx
                intLinY = reEy*reEy + imEy*imEy
                _stokes.arS[i] = intLinX + intLinY
                #_stokes.arS[i + nTot] = intLinX - intLinY
                if(_n_stokes_comp > 1): _stokes.arS[i + nTot] = intLinX - intLinY #OC04052018
                #_stokes.arS[i + nTot2] = -2*(reEx*reEy + imEx*imEy) #check sign
                if(_n_stokes_comp > 2): _stokes.arS[i + nTot2] = 2*(reEx*reEy + imEx*imEy) #OC04052018 #check sign (in SRW for Igor: -2*(ReEX*ReEZ + ImEX*ImEZ))
                #_stokes.arS[i + nTot3] = 2*(-reEx*reEy + imEx*imEy) #check sign
                if(_n_stokes_comp > 3): _stokes.arS[i + nTot3] = 2*(reEx*imEy - imEx*reEy) #OC04052018 #check sign (in SRW for Igor: 2*(-ReEX*ImEZ + ImEX*ReEZ))

                #DEBUG
                #if(isnan(reEx) or isnan(imEx) or isnan(reEy) or isnan(imEy)):
                #    print('reEx=', reEx, ' imEx=', imEx, ' reEy=', reEy, ' imEy=', imEy)
                #END DEBUG
                
            _stokes.mesh.set_from_other(self.mesh)
            #Why this is using self.mesh, whereas if(_stokes.mutual) uses _stokes.mesh? Consider correcting.
            
        else: #calculate Mutual Stokes parameters on the _stokes.mesh
            yNpRes = _stokes.mesh.ny
            yStartRes = _stokes.mesh.yStart
            yStepRes = 0
            if(yNpRes > 1): yStepRes = (_stokes.mesh.yFin - yStartRes)/(yNpRes - 1)

            xNpRes = _stokes.mesh.nx
            xStartRes = _stokes.mesh.xStart
            xStepRes = 0
            if(xNpRes > 1): xStepRes = (_stokes.mesh.xFin - xStartRes)/(xNpRes - 1)

            eNpRes = _stokes.mesh.ne
            eStartRes = _stokes.mesh.eStart
            eStepRes = 0
            if(eNpRes > 1): eStepRes = (_stokes.mesh.eFin - eStartRes)/(eNpRes - 1)

            nTot = eNpRes*xNpRes*yNpRes #OC06052018 (uncommented)
            #nTot = 2*eNpRes*xNpRes*yNpRes #OC04052018 (since <E(r1)E(r2)*> and other Stokes may be complex entities)

            #nTot1 = nTot*nTot
            nTot1 = 2*nTot*nTot #OC06052018
            nTot2 = nTot1*2
            nTot3 = nTot1*3

            yNpWfr = self.mesh.ny
            yStartWfr = self.mesh.yStart
            yStepWfr = 0
            if(yNpWfr > 1): yStepWfr = (self.mesh.yFin - yStartWfr)/(yNpWfr - 1)
            yNpWfr_mi_1 = yNpWfr - 1

            xNpWfr = self.mesh.nx
            xStartWfr = self.mesh.xStart
            xStepWfr = 0
            if(xNpWfr > 1): xStepWfr = (self.mesh.xFin - xStartWfr)/(xNpWfr - 1)
            xNpWfr_mi_1 = xNpWfr - 1

            eNpWfr = self.mesh.ne
            eStartWfr = self.mesh.eStart
            eStepWfr = 0
            if(eNpWfr > 1): eStepWfr = (self.mesh.eFin - eStartWfr)/(eNpWfr - 1)
            eNpWfr_mi_1 = eNpWfr - 1

            perE = 2
            perX = perE*eNpWfr
            perY = perX*xNpWfr

            perXr = perE*eNpRes
            perYr = perX*xNpRes

            nTotAux = nTot*2
            auxArEx = array('f', [0]*nTotAux) #OC06052018
            auxArEy = array('f', [0]*nTotAux)
            #auxArEx = array('f', [0]*nTot) #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)
            #auxArEy = array('f', [0]*nTot)

            #print(perE, perX, perY)

            #ir = 0
            yRes = yStartRes
            for iy in range(yNpRes):
                iyWfr0 = 0
                if(yStepWfr > 0): iyWfr0 = int(trunc((yRes - yStartWfr)/yStepWfr + 1.e-09))
                if((iyWfr0 < 0) or (iyWfr0 > yNpWfr_mi_1)):
                    #_stokes.arS[ir] = 0; _stokes.arS[ir + nTot1] = 0; _stokes.arS[ir + nTot2] = 0; _stokes.arS[ir + nTot3] = 0;
                    #ir += 1;
                    yRes += yStepRes
                    continue
                iyWfr1 = iyWfr0 + 1
                if(iyWfr1 > yNpWfr_mi_1): iyWfr1 = yNpWfr_mi_1
                ty = 0
                if(yStepWfr > 0): ty = (yRes - (yStartWfr + yStepWfr*iyWfr0))/yStepWfr

                iy0_perY = iyWfr0*perY
                iy1_perY = iyWfr1*perY
                iy_perYr = iy*perYr

                #OC21052020
                y_mi_ycE2_d_Ry = 0
                if(_ry_avg != 0):
                    y_mi_yc = yRes - _yc_avg
                    y_mi_ycE2_d_Ry = y_mi_yc*y_mi_yc/_ry_avg

                xRes = xStartRes
                for ix in range(xNpRes):
                    ixWfr0 = 0
                    if(xStepWfr > 0): ixWfr0 = int(trunc((xRes - xStartWfr)/xStepWfr + 1.e-09))
                    if((ixWfr0 < 0) or (ixWfr0 > xNpWfr_mi_1)):
                        #_stokes.arS[ir] = 0; _stokes.arS[ir + nTot1] = 0; _stokes.arS[ir + nTot2] = 0; _stokes.arS[ir + nTot3] = 0;
                        #ir += 1;
                        xRes += xStepRes
                        continue
                    ixWfr1 = ixWfr0 + 1
                    if(ixWfr1 > xNpWfr_mi_1): ixWfr1 = xNpWfr_mi_1
                    tx = 0
                    if(xStepWfr > 0): tx = (xRes - (xStartWfr + xStepWfr*ixWfr0))/xStepWfr

                    ix0_perX = ixWfr0*perX
                    ix1_perX = ixWfr1*perX
                    ix_perXr = ix*perXr

                    #OC21052020
                    x_mi_xcE2_d_Rx = 0
                    if(_rx_avg != 0):
                        x_mi_xc = xRes - _xc_avg
                        x_mi_xcE2_d_Rx = x_mi_xc*x_mi_xc/_rx_avg

                    quadTermBare = x_mi_xcE2_d_Rx + y_mi_ycE2_d_Ry #OC21052020
     
                    eRes = eStartRes
                    for ie in range(eNpRes):
                        ieWfr0 = 0
                        if(eStepWfr > 0): ieWfr0 = int(trunc((eRes - eStartWfr)/eStepWfr + 1.e-09))
                        if((ieWfr0 < 0) or (ieWfr0 > eNpWfr_mi_1)):
                            #_stokes.arS[ir] = 0; _stokes.arS[ir + nTot1] = 0; _stokes.arS[ir + nTot2] = 0; _stokes.arS[ir + nTot3] = 0;
                            #ir += 1;
                            eRes += eStepRes
                            continue
                        ieWfr1 = ieWfr0 + 1
                        if(ieWfr1 > eNpWfr_mi_1): ieWfr1 = eNpWfr_mi_1
                        te = 0
                        if(eStepWfr > 0): te = (eRes - (eStartWfr + eStepWfr*ieWfr0))/eStepWfr

                        ie0_perE = ieWfr0*perE
                        ie1_perE = ieWfr1*perE
                        ie_perE = ie*perE

                        ofstR = ie_perE + ix_perXr + iy_perYr
                        ofstR_p_1 = ofstR + 1 #OC21052020
                                
                        ofst000 = ie0_perE + ix0_perX + iy0_perY
                        ofst100 = ie1_perE + ix0_perX + iy0_perY
                        ofst010 = ie0_perE + ix1_perX + iy0_perY
                        ofst001 = ie0_perE + ix0_perX + iy1_perY
                        ofst110 = ie1_perE + ix1_perX + iy0_perY
                        ofst101 = ie1_perE + ix0_perX + iy1_perY
                        ofst011 = ie0_perE + ix1_perX + iy1_perY
                        ofst111 = ie1_perE + ix1_perX + iy1_perY
                                    
                        a000 = self.arEx[ofst000]#; print(a000)
                        f100 = self.arEx[ofst100]
                        f010 = self.arEx[ofst010]
                        f001 = self.arEx[ofst001]
                        f110 = self.arEx[ofst110]
                        f101 = self.arEx[ofst101]
                        f011 = self.arEx[ofst011]
                        f111 = self.arEx[ofst111]
                        a100 = f100 - a000
                        a010 = f010 - a000
                        a001 = f001 - a000
                        a110 = a000 - f010 - f100 + f110
                        a101 = a000 - f001 - f100 + f101
                        a011 = a000 - f001 - f010 + f011
                        a111 = f001 + f010 - f011 + f100 - f101 - f110 + f111 - a000
                        #auxArEx[ir] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        auxArEx[ofstR] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty

                        a000 = self.arEx[ofst000 + 1]
                        f100 = self.arEx[ofst100 + 1]
                        f010 = self.arEx[ofst010 + 1]
                        f001 = self.arEx[ofst001 + 1]
                        f110 = self.arEx[ofst110 + 1]
                        f101 = self.arEx[ofst101 + 1]
                        f011 = self.arEx[ofst011 + 1]
                        f111 = self.arEx[ofst111 + 1]
                        a100 = f100 - a000
                        a010 = f010 - a000
                        a001 = f001 - a000
                        a110 = a000 - f010 - f100 + f110
                        a101 = a000 - f001 - f100 + f101
                        a011 = a000 - f001 - f010 + f011
                        a111 = f001 + f010 - f011 + f100 - f101 - f110 + f111 - a000
                        #auxArEx[ir + 1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        #auxArEx[ofstR + 1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        auxArEx[ofstR_p_1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty #OC21052020

                        a000 = self.arEy[ofst000]
                        f100 = self.arEy[ofst100]
                        f010 = self.arEy[ofst010]
                        f001 = self.arEy[ofst001]
                        f110 = self.arEy[ofst110]
                        f101 = self.arEy[ofst101]
                        f011 = self.arEy[ofst011]
                        f111 = self.arEy[ofst111]
                        a100 = f100 - a000
                        a010 = f010 - a000
                        a001 = f001 - a000
                        a110 = a000 - f010 - f100 + f110
                        a101 = a000 - f001 - f100 + f101
                        a011 = a000 - f001 - f010 + f011
                        a111 = f001 + f010 - f011 + f100 - f101 - f110 + f111 - a000
                        #auxArEy[ir] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        auxArEy[ofstR] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                                    
                        a000 = self.arEy[ofst000 + 1]
                        f100 = self.arEy[ofst100 + 1]
                        f010 = self.arEy[ofst010 + 1]
                        f001 = self.arEy[ofst001 + 1]
                        f110 = self.arEy[ofst110 + 1]
                        f101 = self.arEy[ofst101 + 1]
                        f011 = self.arEy[ofst011 + 1]
                        f111 = self.arEy[ofst111 + 1]
                        a100 = f100 - a000
                        a010 = f010 - a000
                        a001 = f001 - a000
                        a110 = a000 - f010 - f100 + f110
                        a101 = a000 - f001 - f100 + f101
                        a011 = a000 - f001 - f010 + f011
                        a111 = f001 + f010 - f011 + f100 - f101 - f110 + f111 - a000
                        #auxArEy[ir + 1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        #auxArEy[ofstR + 1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        auxArEy[ofstR_p_1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty #OC21052020

                        #OC21052020
                        if(quadTermBare != 0): #Removing the Quadratic Phase Term if _rx_avg != 0 or _ry_avg != 0
                            #pi_d_lamb = eRes*2.53384080189e+06
                            quadPhTerm = -eRes*2.53384080189e+06*quadTermBare
                            cosPhTerm = cos(quadPhTerm)
                            sinPhTerm = sin(quadPhTerm)
                            reEx = auxArEx[ofstR]
                            imEx = auxArEx[ofstR_p_1]
                            auxArEx[ofstR] = reEx*cosPhTerm - imEx*sinPhTerm
                            auxArEx[ofstR_p_1] = imEx*cosPhTerm + reEx*sinPhTerm
                            reEy = auxArEy[ofstR]
                            imEy = auxArEy[ofstR_p_1]
                            auxArEy[ofstR] = reEy*cosPhTerm - imEy*sinPhTerm
                            auxArEy[ofstR_p_1] = imEy*cosPhTerm + reEy*sinPhTerm

                        #ir += 2
                        eRes += eStepRes
                    xRes += xStepRes
                yRes += yStepRes

            #DEBUG
            #print('SRWLWfr::calc_stokes: eNpRes=', eNpRes, ' xNpRes=', xNpRes, ' yNpRes=', yNpRes)
            #print(' nTot=', nTot, ' nTot1=', nTot1, ' nTot2=', nTot2, ' nTot3=', nTot3)
            #END DEBUG

            perX = perE*eNpRes
            perY = perX*xNpRes
            ir = 0

            for iy in range(yNpRes):
                iy_perY = iy*perY
                for iyp in range(yNpRes):
                    iyp_perY = iyp*perY
                    for ix in range(xNpRes):
                        ix_perX = ix*perX
                        ix_perX_p_iy_perY = ix_perX + iy_perY
                        for ixp in range(xNpRes):
                            ixp_perX = ixp*perX
                            ixp_perX_p_iyp_perY = ixp_perX + iyp_perY
                            for ie in range(eNpRes):
                                ie_perE = ie*perE
                                ie_perE_p_ix_perX_p_iy_perY = ie_perE + ix_perX_p_iy_perY
                                reEx = auxArEx[ie_perE_p_ix_perX_p_iy_perY]
                                imEx = auxArEx[ie_perE_p_ix_perX_p_iy_perY + 1]
                                reEy = auxArEy[ie_perE_p_ix_perX_p_iy_perY]
                                imEy = auxArEy[ie_perE_p_ix_perX_p_iy_perY + 1]
                                for iep in range(eNpRes):
                                    iep_perE = iep*perE
                                    iep_perE_p_ixp_perX_p_iyp_perY = iep_perE + ixp_perX_p_iyp_perY
                                    reExT = auxArEx[iep_perE_p_ixp_perX_p_iyp_perY]
                                    imExT = auxArEx[iep_perE_p_ixp_perX_p_iyp_perY + 1]
                                    reEyT = auxArEy[iep_perE_p_ixp_perX_p_iyp_perY]
                                    imEyT = auxArEy[iep_perE_p_ixp_perX_p_iyp_perY + 1]

                                    #OC04052018 (commented-out)
                                    #intLinX = reEx*reExT + imEx*imExT
                                    #intLinY = reEy*reEyT + imEy*imEyT #; print(intLinX, intLinY)
                                    #_stokes.arS[ir] = intLinX + intLinY
                                    #_stokes.arS[ir + nTot1] = intLinX - intLinY #check sign
                                    #_stokes.arS[ir + nTot2] = -reEx*reEyT - reExT*reEy - imEx*imEyT - imExT*imEy #-2*(reEx*reEy + imEx*imEy) #check sign
                                    #_stokes.arS[ir + nTot3] = -reEx*reEyT - reExT*reEy + imEx*imEyT + imExT*imEy #2*(-reEx*reEy + imEx*imEy) #check sign

                                    #OC04052018
                                    reMI_s0_X = reEx*reExT + imEx*imExT
                                    reMI_s0_Y = reEy*reEyT + imEy*imEyT
                                    imMI_s0_X = imEx*reExT - reEx*imExT
                                    imMI_s0_Y = imEy*reEyT - reEy*imEyT
                                    _stokes.arS[ir] = reMI_s0_X + reMI_s0_Y #Re(s0)
                                    ir_p_1 = ir + 1
                                    _stokes.arS[ir_p_1] = imMI_s0_X + imMI_s0_Y #Im(s0)
                                    if(_n_stokes_comp > 1):
                                        _stokes.arS[ir + nTot1] = reMI_s0_X - reMI_s0_Y #Re(s1), see MutualStokes.nb and LectureNotes notebook
                                        _stokes.arS[ir_p_1 + nTot1] = imMI_s0_X - imMI_s0_Y #Im(s1)
                                    if(_n_stokes_comp > 2):
                                        #DEBUG
                                        #print(' ir + nTot2=', ir + nTot2)
                                        #END DEBUG
                                        _stokes.arS[ir + nTot2] = imExT*imEy + imEx*imEyT + reExT*reEy + reEx*reEyT #Re(s2)
                                        _stokes.arS[ir_p_1 + nTot2] = imEy*reExT - imEyT*reEx - imExT*reEy + imEx*reEyT #Im(s2)
                                    if(_n_stokes_comp > 3):
                                        _stokes.arS[ir + nTot3] = imEyT*reEx + imEy*reExT - imExT*reEy - imEx*reEyT #Re(s3)
                                        _stokes.arS[ir_p_1 + nTot3] = imEx*imEyT - imExT*imEy - reExT*reEy + reEx*reEyT #Im(s3)

                                    #ir += 1
                                    ir += 2 #OC04052018 (since <E(r1)E(r2)*> and other Stokes are complex entities)
            del auxArEx
            del auxArEy
        
    def sim_src_offset(self, _dx, _dxp, _dy, _dyp, _move_mesh=False, _copy=False): #OC07062022
        '''Simulate offset of the wavefront source in horizontal / vertical position and/or angle (version using opt. elem. (kick and shift))'''

        opAng = SRWLOptAng(_ang_x=_dxp, _ang_y=_dyp)
        pp = [0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]

        wfr = deepcopy(self) if _copy else self

        R = 0.5*(wfr.Rx + wfr.Ry)

        #DEBUG
        #print('sim_src_offset: R = ', R, 'm')
        #END DEBUG

        dx = _dx + _dxp*R
        dy = _dy + _dyp*R

        if _move_mesh:

            opCnt = SRWLOptC([opAng], [pp])
            srwl.PropagElecField(wfr, opCnt)

            wfr.mesh.xStart += dx
            wfr.mesh.xFin += dx
            wfr.mesh.yStart += dy
            wfr.mesh.yFin += dy
            wfr.xc += _dx
            wfr.yc += _dy

        else:
    
            opShift = SRWLOptShift(_shift_x=dx, _shift_y=dy)
            opCnt = SRWLOptC([opShift, opAng], [pp, pp])
            srwl.PropagElecField(wfr, opCnt)

            wfr.xc -= _dxp*R
            wfr.yc -= _dyp*R

        #Updating stat. moments of radiation beam:
        wfr.arMomX[0] += dx
        wfr.arMomX[1] += _dxp
        wfr.arMomX[2] += dy
        wfr.arMomX[3] += _dyp
        wfr.arMomY[0] += dx
        wfr.arMomY[1] += _dxp
        wfr.arMomY[2] += dy
        wfr.arMomY[3] += _dyp

        #Updating stat. momemnts of source:
        wfr.partBeam.partStatMom1.x += _dx 
        wfr.partBeam.partStatMom1.xp += _dxp
        wfr.partBeam.partStatMom1.y += _dy
        wfr.partBeam.partStatMom1.yp += _dyp

        return wfr

#****************************************************************************
class SRWLOpt(object):
    """Optical Element (base class)"""

    def get_orient(self, _e=0): #OC17112019 #To be overridden in derived classes
        tv = [1,0,0]; sv = [0,1,0]; nv = [0,0,1]
        ex = [1,0,0]; ey = [0,1,0]; ez = [0,0,1]
        return [[tv, sv, nv], [ex, ey, ez], [ex, ey, ez]] #[2] is transposed of [ex, ey, ez], to be used for fast space transformation calc.

    def set_rand_par(self, _rand_par): #OC23042020
        """Sets list of params to be eventually randomized in some types of calculations
        :param _rand_par: list of params to be randomized; each element of this list should be: ['param_name', val_avg, val_range, meth]
        """
        #Add checking / parsing _rand_par content here?
        self.RandParam = _rand_par

    def randomize(self): #OC23042020
        """Randomizes parameters of optical element according to self.RandParam to simulate e.g. impact of vibrations on coherence (in P-C calculations)
        """
        if(not hasattr(self, 'RandParam')): return
        if(self.RandParam is None): return
        nPar = len(self.RandParam)
        if(nPar == 0): return
        if(not isinstance(self.RandParam, list)): return

        randMeth = 'uni' #uniform distribution by default

        for i in range(nPar):
            curParData = self.RandParam[i]
            if(not isinstance(curParData, list)): continue
            lenCurParData = len(curParData)
            if(lenCurParData < 2): continue
            curParName = curParData[0]
            if(not hasattr(self, curParName)): continue

            newVal = curParData[1]
            if(lenCurParData > 2):
                randRange = curParData[2]
                if(lenCurParData > 3):
                    randMeth = curParData[3]
                    randMeth = randMeth.lower()
                if(randMeth == 'uni'):
                    randAmp = 0.5*randRange
                    newVal += random.uniform(-randAmp, randAmp)
                elif(randMeth == 'gsn'):
                    randRMS = randRange/2.354820045
                    newVal += random.gauss(0, randRMS)
            setattr(self, curParName, newVal)

class SRWLOptD(SRWLOpt):
    """Optical Element: Drift Space"""
    
    def __init__(self, _L=0, _treat=0):
        """
        :param _L: Length [m]
        :param _treat: switch specifying whether the absolute optical path should be taken into account in radiation phase (=1) or not (=0, default)
        """
        self.L = _L
        self.treat = _treat

class SRWLOptA(SRWLOpt):
    """Optical Element: Aperture / Obstacle"""
    
    def __init__(self, _shape='r', _ap_or_ob='a', _Dx=0, _Dy=0, _x=0, _y=0):
        """
        :param _shape: 'r' for rectangular, 'c' for circular
        :param _ap_or_ob: 'a' for aperture, 'o' for obstacle
        :param _Dx: horizontal transverse dimension [m]; in case of circular aperture, only Dx is used for diameter
        :param _Dy: vertical transverse dimension [m]; in case of circular aperture, Dy is ignored
        :param _x: horizontal transverse coordinate of center [m]
        :param _y: vertical transverse coordinate of center [m]
        """
        self.shape = _shape #'r' for rectangular, 'c' for circular
        self.ap_or_ob = _ap_or_ob #'a' for aperture, 'o' for obstacle
        self.Dx = _Dx #transverse dimensions [m]; in case of circular aperture, only Dx is used for diameter
        self.Dy = _Dy
        self.x = _x #transverse coordinates of center [m]
        self.y = _y

class SRWLOptL(SRWLOpt):
    """Optical Element: Thin Lens"""
    
    def __init__(self, _Fx=1e+23, _Fy=1e+23, _x=0, _y=0):
        """
        :param _Fx: focal length in horizontal plane [m]
        :param _Fy: focal length in vertical plane [m]
        :param _x: horizontal coordinate of center [m]
        :param _y: vertical coordinate of center [m]
        """
        self.Fx = _Fx #focal lengths [m]
        self.Fy = _Fy
        self.x = _x #transverse coordinates of center [m]
        self.y = _y

class SRWLOptAng(SRWLOpt):
    """Optical Element: Angle"""
    
    def __init__(self, _ang_x=0, _ang_y=0):
        """
        :param _ang_x: horizontal angle [rad]
        :param _ang_y: vertical angle [rad]
        """
        self.AngX = _ang_x
        self.AngY = _ang_y

class SRWLOptShift(SRWLOpt):
    """Optical Element: Shirt"""
    
    def __init__(self, _shift_x=0, _shift_y=0):
        """
        :param _shift_x: horizontal shift [m]
        :param _shift_y: vertical shift [m]
        """
        self.ShiftX = _shift_x
        self.ShiftY = _shift_y

class SRWLOptZP(SRWLOpt):
    """Optical Element: Thin Lens"""
    
    def __init__(self, _nZones=100, _rn=0.1e-03, _thick=10e-06, _delta1=1e-06, _atLen1=0.1, _delta2=0, _atLen2=1e-06, _x=0, _y=0, _e=0): #OC22062019
    #def __init__(self, _nZones=100, _rn=0.1e-03, _thick=10e-06, _delta1=1e-06, _atLen1=0.1, _delta2=0, _atLen2=1e-06, _x=0, _y=0):
        """
        :param _nZones: total number of zones
        :param _rn: auter zone radius [m]
        :param _thick: thickness [m]
        :param _delta1: refractuve index decrement of the "main" material
        :param _atLen1: attenuation length [m] of the "main" material
        :param _delta2: refractuve index decrement of the "complementary" material
        :param _atLen2: attenuation length [m] of the "complementary" material
        :param _x: horizontal transverse coordinate of center [m]
        :param _y: vertical transverse coordinates of center [m]
        :param _e: average photon energy [eV], active if > 0, to be used for calculation corrections to zone radii
        """
        self.nZones = _nZones #total number of zones
        self.rn = _rn #auter zone radius [m]
        self.thick = _thick #thickness [m]
        self.delta1 = _delta1 #refractuve index decrement of the "main" material
        self.delta2 = _delta2 #refractuve index decrement of the "complementary" material
        self.atLen1 = _atLen1 #attenuation length [m] of the "main" material
        self.atLen2 = _atLen2 #attenuation length [m] of the "complementary" material
        self.x = _x #transverse coordinates of center [m]
        self.y = _y
        self.e0 = _e #average photon energy [eV], to be used for calculation corrections to zone radii

    def get_F(self, _e): #OC22062019
        """Estimate focal length"""

        mult = _Light_eV_mu*1.e-06 #photon energy <-> wavelength
        lamb = mult/_e

        aux = (self.rn)*(self.rn)*(1 - 1/self.nZones)
        if((hasattr(self, e0))):
           if(self.e0 > 0):
               lamb0 = mult/self.e0
               aux -= 0.25*lamb0*lamb0*(self.nZones - 1)
        rn_mi_1 = sqrt(aux)
        two_drn = self.rn - rn_mi_1
        aux = lamb/two_drn
        return (two_drn*self.rn/lamb)*sqrt(1 - aux*aux)

class SRWLOptWG(SRWLOpt):
    """Optical Element: Waveguide"""
    
    def __init__(self, _L=1, _Dx=10e-03, _Dy=10e-03, _x=0, _y=0):
        """
        :param _L: length [m]
        :param _Dx: horizontal transverse dimension [m]
        :param _Dy: vertical transverse dimension [m]
        :param _x: horizontal transverse coordinate of center [m]
        :param _y: vertical transverse coordinate of center [m]
        """
        self.L = _L #length [m]
        self.Dx = _Dx #transverse dimensions [m]
        self.Dy = _Dy
        self.x = _x #transverse coordinates of center [m]
        self.y = _y

class SRWLOptT(SRWLOpt):
    """Optical Element: Transmission (generic)"""
    
    def __init__(self, _nx=1, _ny=1, _rx=1e-03, _ry=1e-03, _arTr=None, _extTr=0, _Fx=1e+23, _Fy=1e+23, _x=0, _y=0, _ne=1, _eStart=0, _eFin=0, _alloc_base=[0]): #OC14042019
    #def __init__(self, _nx=1, _ny=1, _rx=1e-03, _ry=1e-03, _arTr=None, _extTr=0, _Fx=1e+23, _Fy=1e+23, _x=0, _y=0, _ne=1, _eStart=0, _eFin=0):
        """
        :param _nx: number of transmission data points in the horizontal direction
        :param _ny: number of transmission data points in the vertical direction
        :param _rx: range of the horizontal coordinate [m] for which the transmission is defined
        :param _ry: range of the vertical coordinate [m] for which the transmission is defined
        :param _arTr: complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse coordinates
        :param _extTr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1)
        :param _Fx: estimated focal length in the horizontal plane [m]
        :param _Fy: estimated focal length in the vertical plane [m]
        :param _x: horizontal transverse coordinate of center [m]
        :param _y: vertical transverse coordinate of center [m]
        :param _ne: number of transmission data points vs photon energy
        :param _eStart: initial value of photon energy
        :param _eFin: final value of photon energy
        """
        
        self.arTr = _arTr #complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse position
        if((_arTr is None) or ((len(_arTr) != _ne*_nx*_ny*2) and (_ne*_nx*_ny > 0))):
            self.allocate(_ne, _nx, _ny, _alloc_base) #OC14042019
            #self.allocate(_ne, _nx, _ny)
            #print(_ne, _nx, _ny)

        #self.ne = _ne #number of transmission data points vs photon energy
        #self.nx = _nx #numbers of transmission data points in the horizontal and vertical directions
        #self.ny = _ny
        #self.eStart = _eStart #initial and final values of photon energy
        #self.eFin = _eFin
        #self.rx = _rx #ranges of horizontal and vertical coordinates [m] for which the transmission is defined
        #self.ry = _ry

        halfRangeX = 0.5*_rx
        halfRangeY = 0.5*_ry

        if(not hasattr(self, 'mesh')): #OC10112018
            self.mesh = SRWLRadMesh(_eStart, _eFin, _ne, _x - halfRangeX, _x + halfRangeX, _nx, _y - halfRangeY, _y + halfRangeY, _ny)
        else:
            self.mesh.eStart = _eStart
            self.mesh.eFin = _eFin
            self.mesh.xStart = _x - halfRangeX
            self.mesh.xFin = _x + halfRangeX
            self.mesh.yStart = _y - halfRangeY
            self.mesh.yFin = _y + halfRangeY

        self.extTr = _extTr #0- transmission outside the grid/mesh is zero; 1- it is same as on boundary
        self.Fx = _Fx #estimated focal lengths [m]
        self.Fy = _Fy
        
        #self.x = _x #transverse coordinates of center [m]
        #self.y = _y
        #if _ne > 1: _Fx, _Fy should be arrays vs photon energy?

    def allocate(self, _ne, _nx, _ny, _alloc_base=[0]): #OC14042019
    #def allocate(self, _ne, _nx, _ny):
        #self.ne = _ne
        #self.nx = _nx
        #self.ny = _ny

        if(hasattr(self, 'mesh')):
            self.mesh.ne = _ne
            self.mesh.nx = _nx
            self.mesh.ny = _ny            
        else:
            self.mesh = SRWLRadMesh(0, 0, _ne, 0, 0, _nx, 0, 0, _ny)

        nTot = 2*_ne*_nx*_ny #total array length to store amplitude transmission and optical path difference
        #self.arTr = array('d', [0]*nTot)
        
        lenBase = len(_alloc_base) #OC14042019
        if(lenBase > 1): nTot = int(round(nTot/lenBase)) #OC14042019
        self.arTr = srwl_uti_array_alloc('d', nTot, _alloc_base) #OC14042019

    def get_data(self, _typ, _dep=3, _e=0, _x=0, _y=0):
        """Returns Transmission Data Characteristic
        :param _typ: type of transmission characteristic to extract: 1- amplitude transmission, 2- intensity transmission, 3- optical path difference
        :param _dep: type of dependence to extract: 0- vs photon energy, 1- vs horizontal position, 2- vs vertical position, 3- vs hor. & vert. positions
        :param _e: photon energy [eV] (to keep fixed)
        :param _x: horizontal position [m] (to keep fixed)
        :param _y: vertical position [m] (to keep fixed)
        """
        nTot = self.mesh.ne*self.mesh.nx*self.mesh.ny
        arAux = array('d', [0]*nTot)
        for i in range(nTot): #put all data into one column using "C-alignment" as a "flat" 1D array
            tr = 0
            if((_typ == 1) or (_typ == 2)): #amplitude or intensity transmission
                tr = self.arTr[i*2]
                if(_typ == 2): #intensity transmission
                    tr *= tr
            else: #optical path difference
                tr = self.arTr[i*2 + 1]
            arAux[i] = tr
        if (_dep == 3) and (self.mesh.ne == 1): return arAux
        #print('total extract passed')
        
        arOut = None
        xStep = 0
        if self.mesh.nx > 1: xStep = (self.mesh.xFin - self.mesh.xStart)/(self.mesh.nx - 1)
        yStep = 0
        if self.mesh.ny > 1: yStep = (self.mesh.yFin - self.mesh.yStart)/(self.mesh.ny - 1)
        inperpOrd = 1 #inperpolation order, up to 3
        if _dep == 0: #dependence vs photon energy
            arOut = array('d', [0]*self.mesh.ne)
            for ie in range(self.mesh.ne):
                #arOut[ie] = srwl_uti_interp_2d(_x, _y, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                arOut[ie] = uti_math.interp_2d(_x, _y, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
        else:
            ie = 0
            if self.mesh.ne > 1:
                if _e >= self.mesh.eFin: ie = self.mesh.ne - 1
                elif _e > self.mesh.eStart:
                    eStep = (self.mesh.eFin - self.mesh.eStart)/(self.mesh.ne - 1)
                    ie = int(round((_e - self.mesh.eStart)/eStep))
            #print(ie)
            if _dep == 1: #dependence vs horizontal position
                arOut = array('d', [0]*self.mesh.nx)
                xx = self.mesh.xStart
                for ix in range(self.mesh.nx):
                    #arOut[ix] = srwl_uti_interp_2d(xx, _y, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                    arOut[ix] = uti_math.interp_2d(xx, _y, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                    xx += xStep
            elif _dep == 2: #dependence vs vertical position
                arOut = array('d', [0]*self.mesh.ny)
                yy = self.mesh.yStart
                for iy in range(self.mesh.ny):
                    #arOut[iy] = srwl_uti_interp_2d(_x, yy, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                    arOut[iy] = uti_math.interp_2d(_x, yy, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                    yy += yStep
            elif _dep == 3: #dependence vs horizontal and vertical position
                nTot = self.mesh.nx*self.mesh.ny
                arOut = array('d', [0]*nTot)
                yy = self.mesh.yStart
                i = 0
                for iy in range(self.mesh.ny):
                    xx = self.mesh.xStart
                    for ix in range(self.mesh.nx):
                        #arOut[i] = srwl_uti_interp_2d(xx, yy, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                        arOut[i] = uti_math.interp_2d(xx, yy, self.mesh.xStart, xStep, self.mesh.nx, self.mesh.yStart, yStep, self.mesh.ny, arAux, inperpOrd, self.mesh.ne, ie)
                        i += 1
                        xx += xStep
                    yy += yStep
        del arAux
        #print(len(arOut))
        return arOut

    def randomize(self): #OC28042020 (overwriting base-class function)
        """Randomizes parameters of optical element according to self.RandParam to simulate e.g. impact of vibrations on coherence (in P-C calculations)
        """
        if(not hasattr(self, 'RandParam')): return
        if(self.RandParam is None): return
        nPar = len(self.RandParam)
        if(nPar == 0): return
        if(not isinstance(self.RandParam, list)): return

        randMeth = 'uni' #uniform distribution by default

        for i in range(nPar):
            curParData = self.RandParam[i]
            if(not isinstance(curParData, list)): continue
            lenCurParData = len(curParData)
            if(lenCurParData < 2): continue
            curParName = curParData[0]

            if((curParName == 'x') or (curParName == 'y')):
            #if(not hasattr(self, curParName)): continue
                newVal = curParData[1]
                if(lenCurParData > 2):
                    randRange = curParData[2]
                    if(lenCurParData > 3):
                        randMeth = curParData[3]
                        randMeth = randMeth.lower()
                    if(randMeth == 'uni'):
                        randAmp = 0.5*randRange
                        newVal += random.uniform(-randAmp, randAmp)
                    elif(randMeth == 'gsn'):
                        randRMS = randRange/2.354820045
                        newVal += random.gauss(0, randRMS)

                #setattr(self, curParName, newVal)
                if(curParName == 'x'):
                    xHalfRange = 0.5*(self.mesh.xFin - self.mesh.xStart)
                    self.mesh.xStart = newVal - xHalfRange
                    self.mesh.xFin = newVal + xHalfRange
                elif(curParName == 'y'):
                    yHalfRange = 0.5*(self.mesh.yFin - self.mesh.yStart)
                    self.mesh.yStart = newVal - yHalfRange
                    self.mesh.yFin = newVal + yHalfRange
                elif(curParName == 'e'):
                    eHalfRange = 0.5*(self.mesh.eFin - self.mesh.eStart)
                    self.mesh.eStart = newVal - eHalfRange
                    self.mesh.eFin = newVal + eHalfRange

class SRWLOptMir(SRWLOpt):
    """Optical Element: Mirror (focusing)"""

    def set_dim_sim_meth(self, _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=500, _nps=500, _treat_in_out=1, _ext_in=0, _ext_out=0):
        """Sets Mirror Dimensions, Aperture Shape and its simulation method
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _npt: number of mesh points to represent mirror in tangential direction (used for "thin" approximation)
        :param _nps: number of mesh points to represent mirror in sagital direction (used for "thin" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after this the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        """
        if((_sim_meth < 1) or (_sim_meth > 2)):
            raise Exception("Simulation method is not specified correctly (should be 1 for \"thin\", 2 for \"thick\" element approximation)")
        self.dt = _size_tang
        self.ds = _size_sag
        self.apShape = _ap_shape
        self.meth = _sim_meth
        self.npt = _npt
        self.nps = _nps
        self.treatInOut = _treat_in_out
        self.extIn = _ext_in
        self.extOut = _ext_out
        self.Fx = 0 #i.e. focal lengthes are not set
        self.Fy = 0

    def set_reflect(self, _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=0, _ph_en_fin=0, _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """Sets Mirror Reflectivity
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat array complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)
        """
        nTot = int(_n_ph_en*_n_ang*_n_comp*2)
        if(nTot < 2):
            raise Exception("Incorrect Reflectivity array parameters")
        _n_comp = int(_n_comp)
        if((_n_comp < 1) or (_n_comp > 2)):
            raise Exception("Number of reflectivity coefficient components can be 1 or 2")

        self.arRefl = None #OC12082018
        if((_refl is not None) and (_refl != 1)): #OC12082018 
            if(not(isinstance(_refl, list) or isinstance(_refl, array))):
                self.arRefl = array('d', [_refl]*nTot)
                for i in range(int(round(nTot/2))):
                    i2 = i*2
                    self.arRefl[i2] = _refl
                    self.arRefl[i2 + 1] = 0
            else:
                self.arRefl = _refl

        #DEBUG
        #print('Reflection:', _refl, self.arRefl)
        
        self.reflNumPhEn = int(_n_ph_en)
        self.reflNumAng = int(_n_ang)
        self.reflNumComp = _n_comp
        self.reflPhEnStart = _ph_en_start
        self.reflPhEnFin = _ph_en_fin
        self.reflPhEnScaleType = _ph_en_scale_type
        self.reflAngStart = _ang_start
        self.reflAngFin = _ang_fin
        self.reflAngScaleType = _ang_scale_type

    def set_orient(self, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0):
        """Defines Mirror Orientation in the frame of the incident photon beam
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        """
        self.nvx = _nvx
        self.nvy = _nvy
        self.nvz = _nvz
        self.tvx = _tvx
        self.tvy = _tvy
        self.x = _x
        self.y = _y
        self.Fx = 0 #i.e. focal lengths are (re-)set, because changing orientation affects them
        self.Fy = 0

    def set_all(self, _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=100, _nps=100, _treat_in_out=1, _ext_in=0, _ext_out=0,
                _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after that the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)      
        """

        self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

    def get_orient(self, _e=0): #OC18112019

        nv = [self.nvx, self.nvy, self.nvz]
        tv = [self.tvx, self.tvy, -(self.nvx*self.tvx + self.nvy*self.tvy)/self.nvz]; 
        sv = uti_math.vect3_prod_v(nv, tv)

        #Base vectors of the output beam (in the frame of incident beam)
        ezIn = [0,0,1]
        mult = -2*uti_math.vect_prod_s(ezIn, nv)
        ez = [(ezIn[i] + mult*nv[i]) for i in range(3)]
        ey = None
        ex = None

        ezIn_ez_scal = uti_math.vect_prod_s(ezIn, ez)
        angTol = 1e-07
        if(abs(ezIn_ez_scal + 1.) <= angTol):
            ex = [-1,0,0]
            ey = [0,1,0]
        else:
            ezIn_ez_vect = uti_math.vect3_prod_v(ezIn, ez)
            vAxRot = uti_math.vect_normalize(copy(ezIn_ez_vect))
            sinAng = uti_math.vect_prod_s(ezIn_ez_vect, vAxRot)
            cosAng = ezIn_ez_scal
            angRot = 0 
            #The following will determine -pi < angRot <= pi
            if(cosAng < 0):
                if(sinAng < 0): angRot = -pi + asin(-sinAng)
                else: angRot = pi - asin(sinAng)
            else: angRot = asin(sinAng)

            Mrot = uti_math.trf_rotation(vAxRot, angRot, [0]*3)[0] #Matrix of Rotation
            ex = uti_math.matr_prod(Mrot, [1,0,0])
            ey = uti_math.matr_prod(Mrot, [0,1,0])
            
            #DEBUG
            #print('Mirror Loc. Frame nv=', nv, ' ez=', ez)
            #print('norm_ex=', uti_math.vect_norm(ex), 'norm_ey=', uti_math.vect_norm(ey), 'norm_ez=', uti_math.vect_norm(ez))
            #print('Mirror Mrot and [ex, ey, ez]:')
            #print(Mrot)
            #print([ex, ey, ez])

        return [[tv, sv, nv], [ex, ey, ez], Mrot] #[2] is transposed of [ex, ey, ez], to be used for fast space transformation calc.

class SRWLOptMirPl(SRWLOptMir):
    """Optical Element: Mirror: Plane"""
    
    def __init__(self, 
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=100, _nps=100, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after that the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)      
        """
        #There are no other members, except for those of the base class.
        #Finishing of the mirror setup requires calling these 3 functions (with their required arguments):
        #self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        #self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        #self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)
        self.set_all(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out,
                     _nvx, _nvy, _nvz, _tvx, _tvy, _x, _y,
                     _refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptMirEl(SRWLOptMir):
    """Optical Element: Mirror: Elliptical
       NOTE: in the Local frame of the Mirror tangential direction is X, saggital Y, mirror normal is along Z"""
    
    def __init__(self, _p=1, _q=1, _ang_graz=1e-03, _r_sag=1.e+23,
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=500, _nps=500, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _p: distance from first focus (\"source\") to mirror center [m]
        :param _q: distance from mirror center to second focus (\"image\") [m]
        :param _ang_graz: grazing angle at mirror center at perfect orientation [rad]
        :param _r_sag: sagital radius of curvature at mirror center [m]
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _npt: number of mesh points to represent mirror in tangential direction (used for "thin" approximation)
        :param _nps: number of mesh points to represent mirror in sagital direction (used for "thin" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after this the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)      
        """

        self.p = _p
        self.q = _q
        self.angGraz = _ang_graz
        self.radSag = _r_sag
        
        #finishing of the mirror setup requires calling these 3 functions (with their required arguments):
        #self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        #self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        #self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)
        self.set_all(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out,
                     _nvx, _nvy, _nvz, _tvx, _tvy, _x, _y,
                     _refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptMirHyp(SRWLOptMir):
    """Optical Element: Mirror: Hyperboloid #TW24012024
       NOTE: in the Local frame of the Mirror tangential direction is X, saggital Y, mirror normal is along Z"""
    
    def __init__(self, _p=1, _q=1, _ang_graz=1e-03, _r_sag=1.e+23,
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=500, _nps=500, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _p: distance from first focus (\"source\") to mirror center [m]
        :param _q: distance from mirror center to second focus (\"image\") [m]
        :param _ang_graz: grazing angle at mirror center at perfect orientation [rad]
        :param _r_sag: sagital radius of curvature at mirror center [m]
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _npt: number of mesh points to represent mirror in tangential direction (used for "thin" approximation)
        :param _nps: number of mesh points to represent mirror in sagital direction (used for "thin" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after this the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)      
        """

        self.p = _p
        self.q = _q
        self.angGraz = _ang_graz
        self.radSag = _r_sag
        
        #finishing of the mirror setup requires calling these 3 functions (with their required arguments):
        #self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        #self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        #self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)
        self.set_all(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out,
                     _nvx, _nvy, _nvz, _tvx, _tvy, _x, _y,
                     _refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptMirPar(SRWLOptMir):
    """Optical Element: Mirror: Paraboloid
       NOTE: in the Local frame of the Mirror tangential direction is X, saggital Y, mirror normal is along Z"""
    
    def __init__(self, _f=1, _uc='f', _ang_graz=1e-03, _r_sag=1.e+23,
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=500, _nps=500, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _f: focal length [m]
        :param _uc: use case: if _uc == 'f': assumes focusing mirror with source at infinity; if _uc == 'c': assumes collimating mirror with image at infinity
        :param _ang_graz: grazing angle at mirror center at perfect orientation [rad]
        :param _r_sag: sagital radius of curvature at mirror center [m]
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _npt: number of mesh points to represent mirror in tangential direction (used for "thin" approximation)
        :param _nps: number of mesh points to represent mirror in sagital direction (used for "thin" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after this the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)      
        """

        self.f = _f
        self.uc = _uc
        self.angGraz = _ang_graz
        self.radSag = _r_sag
        
        #finishing of the mirror setup requires calling these 3 functions (with their required arguments):
        #self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        #self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        #self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)
        self.set_all(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out,
                     _nvx, _nvy, _nvz, _tvx, _tvy, _x, _y,
                     _refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptMirSph(SRWLOptMir):
    """Optical Element: Mirror: Spherical"""

    def __init__(self, _r=1.,
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=500, _nps=500, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _r: radius of surface curvature [m]
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _npt: number of mesh points to represent mirror in tangential direction (used for "thin" approximation)
        :param _nps: number of mesh points to represent mirror in sagital direction (used for "thin" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after this the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)              
        """
        
        self.rad = _r
        
        #finishing of the mirror setup requires calling these 3 functions (with their required arguments):
        self.set_all(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out,
                     _nvx, _nvy, _nvz, _tvx, _tvy, _x, _y,
                     _refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptMirTor(SRWLOptMir):
    """Optical Element: Mirror: Toroid
       NOTE: in the Local frame of the Mirror tangential direction is X, saggital Y, mirror normal is along Z"""
    
    def __init__(self, _rt=1, _rs=1,
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=500, _nps=500, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        
        """
        :param _rt: tangential (major) radius [m]
        :param _rs: sagittal (minor) radius [m]
        :param _size_tang: size in tangential direction [m]
        :param _size_sag: size in sagital direction [m]
        :param _ap_shape: shape of aperture in local frame ('r' for rectangular, 'e' for elliptical)
        :param _sim_meth: simulation method (1 for "thin" approximation, 2 for "thick" approximation)
        :param _npt: number of mesh points to represent mirror in tangential direction (used for "thin" approximation)
        :param _nps: number of mesh points to represent mirror in sagital direction (used for "thin" approximation)
        :param _treat_in_out: switch specifying how to treat input and output wavefront before and after the main propagation through the optical element:
                0- assume that the input wavefront is defined in the plane before the optical element, and the output wavefront is required in a plane just after the element;
                1- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center;
                2- assume that the input wavefront is defined in the plane at the optical element center and the output wavefront is also required at the element center; however, before the propagation though the optical element, the wavefront should be propagated through a drift back to a plane just before the optical element, then a special propagator will bring the wavefront to a plane at the optical element exit, and after this the wavefront will be propagated through a drift back to the element center;
        :param _ext_in: optical element extent on the input side, i.e. distance between the input plane and the optical center (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters
        :param _ext_out: optical element extent on the output side, i.e. distance between the optical center and the output plane (positive, in [m]) to be used at wavefront propagation manipulations; if 0, this extent will be calculated internally from optical element parameters        
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        :param _refl: reflectivity coefficient to set (can be one number or C-aligned flat complex array vs photon energy vs grazing angle vs component (sigma, pi))
        :param _n_ph_en: number of photon energy values for which the reflectivity coefficient is specified
        :param _n_ang: number of grazing angle values for which the reflectivity coefficient is specified
        :param _n_comp: number of electric field components for which the reflectivity coefficient is specified (can be 1 or 2)
        :param _ph_en_start: initial photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_fin: final photon energy value for which the reflectivity coefficient is specified
        :param _ph_en_scale_type: photon energy sampling type ('lin' for linear, 'log' for logarithmic)
        :param _ang_start: initial grazing angle value for which the reflectivity coefficient is specified
        :param _ang_fin: final grazing angle value for which the reflectivity coefficient is specified
        :param _ang_scale_type: angle sampling type ('lin' for linear, 'log' for logarithmic)                      
        """
        
        self.radTan = _rt
        self.radSag = _rs
        
        #finishing of the mirror setup requires calling these 3 functions (with their required arguments):
        self.set_all(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out,
                     _nvx, _nvy, _nvz, _tvx, _tvy, _x, _y,
                     _refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptG(SRWLOpt):
    """Optical Element: Grating"""
    
    #def __init__(self, _mirSub, _m=1, _grDen=100, _grDen1=0, _grDen2=0, _grDen3=0, _grDen4=0, _grAng=0):
    def __init__(self, _mirSub, _m=1, _grDen=100, _grDen1=0, _grDen2=0, _grDen3=0, _grDen4=0, _grAng=0, _e_avg=0, _cff=None, _ang_graz=0, _ang_roll=0): #OC20112019
        """
        :param _mirSub: SRWLOptMir (or derived) type object defining substrate of the grating
        :param _m: output (diffraction) order
        :param _grDen: groove density [lines/mm] (coefficient a0 in the polynomial groove density: a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4)
        :param _grDen1: groove density polynomial coefficient a1 [lines/mm^2]
        :param _grDen2: groove density polynomial coefficient a2 [lines/mm^3]
        :param _grDen3: groove density polynomial coefficient a3 [lines/mm^4]
        :param _grDen4: groove density polynomial coefficient a4 [lines/mm^5]
        :param _grAng: angle between the grove direction and the saggital direction of the substrate [rad] (by default, groves are made along saggital direction (_grAng=0))
        :param _e_avg: average photon energy [eV] the grating should be aligned for (is taken into account if > 0)
        :param _cff: PGM cff parameter, i.e. cos(beta)/cos(alpha); it will be taken into account if _e_avg != 0
        :param _ang_graz: grazing incidence angle [rad] the grating should be aligned for; it will be taken into account if _e_avg != 0 and _cff is None
        :param _ang_roll: roll angle [rad] (i.e. angle of diffraction plane angle rotation about incident beam axis) the grating should be alligned for: it is taken into account when _e_avg != 0; _ang_roll = 0 corresponds to the vertical beam deflection; pi/2 to the horizontal deflection; any value in between is allowed
        """

        if(isinstance(_mirSub, SRWLOptMir) == False):
            raise Exception("Incorrect substrate data submitted to Grating constructor(SRWLOptMir type object is expected for the substrate)") 

        self.mirSub = _mirSub #SRWLOptMir (or derived) type object defining the Grating substrate
        self.m = _m #output order
        self.grDen = _grDen #grove density [lines/mm] (polynomial coefficient a0 in: a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4)
        self.grDen1 = _grDen1 #grove density polynomial coefficient a1 in [lines/mm^2]
        self.grDen2 = _grDen2 #grove density polynomial coefficient a2 in [lines/mm^3]
        self.grDen3 = _grDen3 #grove density polynomial coefficient a3 in [lines/mm^4]
        self.grDen4 = _grDen4 #grove density polynomial coefficient a4 in [lines/mm^5]
        self.grAng = _grAng #angle between the grove direction and the saggital direction of the substrate [rad]
        
        if((_e_avg > 0) and ((_cff is not None) or (_ang_graz != 0))): #OC21112019
            orntData = self.find_orient(_e_avg, _cff, _ang_graz, _ang_roll) #returns [[tv, sv, nv], [ex, ey, ez]]
            opElBaseVects = orntData[0]
            tv = opElBaseVects[0]
            nv = opElBaseVects[2]
            self.set_orient(_nvx=nv[0], _nvy=nv[1], _nvz=nv[2], _tvx=tv[0], _tvy=tv[1])
            #DEBUG
            #print('GRATING: New orientation vector coords:')
            #print('nvx=', nv[0], 'nvy=', nv[1], 'nvz=', nv[2], 'tvx=', tv[0], 'tvy=', tv[1])

    def cff2ang(self, _en, _cff):
        """Calculates Grating grazing angle and deflection angle from photon energy and PGM parameters
        :param _en: radiation photon energy [eV]
        :param _cff: PGM cff parameter, i.e. cos(beta)/cos(alpha)
        :return: grating grazing angle and deflection angle (NOTE: PGM mirror grazing angle aquals to half deflection angle)
        """
        lamb = _Light_eV_mu*1.e-06/_en
        m_lamb_k0 = self.m*lamb*self.grDen*1000
        cff2 = _cff*_cff
        cff2_mi_1 = cff2 - 1
        #OC24032020
        sinBeta = (m_lamb_k0*cff2 - sqrt(m_lamb_k0*m_lamb_k0*cff2 + cff2_mi_1*cff2_mi_1))/cff2_mi_1 if(self.m >= 0) else (m_lamb_k0*cff2 + sqrt(m_lamb_k0*m_lamb_k0*cff2 + cff2_mi_1*cff2_mi_1))/cff2_mi_1 #AH29032020
        #sinBeta = (m_lamb_k0*cff2 - sqrt(m_lamb_k0*m_lamb_k0*cff2 + cff2_mi_1*cff2_mi_1))/cff2_mi_1 if(_m >= 0) else (m_lamb_k0*cff2 + sqrt(m_lamb_k0*m_lamb_k0*cff2 + cff2_mi_1*cff2_mi_1))/cff2_mi_1
        sinAlpha = m_lamb_k0 - sinBeta
        if((sinBeta < -1.) or (sinBeta > 1.) or (sinAlpha < -1.) or (sinAlpha > 1.)): return None, None
        
        alpha = asin(sinAlpha)
        grAng = 0.5*pi - alpha
        beta = asin(sinBeta)
        defAng = pi - alpha + beta
        return grAng, defAng

    def ang2cff(self, _en, _ang_graz): #AH12042019
        """Calculates Grating grazing angle and deflection angle from photon energy and PGM parameters
        :param _en: radiation photon energy [eV]
        :param _ang_graz: grazing incidence angle [rad] the grating should be aligned for; it will be taken into account if _e_avg != 0 and _cff is None
        :return: PGM cff parameter, i.e. cos(beta)/cos(alpha) and deflection angle
        """
        lamb = _Light_eV_mu*1.e-06/_en
        m_lamb_k0 = self.m*lamb*self.grDen*1000
        alpha = 0.5*pi - _ang_graz
        #cff = sqrt((cos(alpha))**2 + 2*m_lamb_k0*sin(alpha) - (m_lamb_k0)**2) / cos(alpha) #OC24032020 (commented-out)
        sinBeta = m_lamb_k0 - sin(alpha)
        if((sinBeta < -1.) or (sinBeta > 1.)): return None, None #OC24032020
        
        beta = asin(sinBeta)
        cff = cos(beta)/cos(alpha) #OC24032020
        
        defAng = pi - alpha + beta
        return cff, defAng

    def angcff2en(self, _cff, _ang_graz): #AH12042019
        """Calculates Grating grazing angle and deflection angle from photon energy and PGM parameters
        :param _cff: PGM cff parameter, i.e. cos(beta)/cos(alpha)_en: radiation photon energy [eV]
        :param _ang_graz: grazing incidence angle [rad] the grating should be aligned for; it will be taken into account if _e_avg != 0 and _cff is None
        :return: radiation photon energy [eV] and deflection angle
        """
        alpha = 0.5*pi - _ang_graz
        #beta = -acos(_cff*cos(alpha))
        cff_cos_alp = _cff*cos(alpha) #OC24032020
        if((cff_cos_alp < -1.) or (cff_cos_alp > 1.)): return None, None
        beta = -acos(cff_cos_alp)
        
        m_lamb_k0 = sin(beta) + sin(alpha)
        lamb = m_lamb_k0 / (self.m*self.grDen*1000)
        en = _Light_eV_mu*1.e-06/lamb
        defAng = pi - alpha + beta
        return en, defAng

    def find_orient(self, _en, _cff=None, _ang_graz=0, _ang_roll=0): #OC21112019
        """Finds optimal crystal orientation in the input beam frame (i.e. surface normal and tangential vectors) and the orientation of the output beam frame (i.e. coordinates of the longitudinal and horizontal vectors in the input beam frame)
        :param _en: photon energy [eV]
        :param _cff: PGM cff parameter, i.e. cos(beta)/cos(alpha); it will be taken into account if _e_avg != 0; if _cff is not None, it dominates over _ang_graz (i.e. it forces its recalculation)
        :param _ang_graz: grazing incidence angle [rad] the grating should be aligned for; it will be taken into account if _e_avg != 0 and _cff is None
        :param _ang_roll: roll angle [rad] (i.e. angle of diffraction plane angle rotation about incident beam axis) the grating should be alligned for: it is taken into account when _e_avg != 0; _ang_roll = 0 corresponds to the vertical beam deflection; pi/2 to the horizontal deflection; any value in between is allowed
        """

        if((_en <= 0) or ((_cff is None) and (_ang_graz == 0))): return None
            
        grAng = _ang_graz
        defAng = 0
        if(_cff is not None): grAng, defAng = self.cff2ang(_en, _cff)
        
        if(grAng <= 0): return None
        
        #DEBUG
        #print('Grating: grazing angle:', grAng)

        alpha = 0.5*pi - grAng
        sinAlpha = sin(alpha)
        cosAlpha = cos(alpha)
        
        nv = [0, sinAlpha, -cosAlpha]
        tv = [0, cosAlpha, sinAlpha] #or it can be [0, -cosAlpha, -sinAlpha] - to re-check assumption used in C++
        
        ezIn = [0,0,1]
        if(_ang_roll != 0):
            MrotZ = uti_math.trf_rotation(ezIn, _ang_roll, [0]*3)[0] #Matrix of Rotation around input beam axis
            nv = uti_math.matr_prod(MrotZ, nv)
            tv = uti_math.matr_prod(MrotZ, tv)
            
        sv = uti_math.vect3_prod_v(nv, tv)
        
        if(defAng == 0):
            lamb = _Light_eV_mu*1.e-06/_en
            m_lamb_k0 = self.m*lamb*self.grDen*1000
            sinBeta = m_lamb_k0 - sinAlpha
            beta = asin(sinBeta)
            defAng = pi - alpha + beta
            
        vAxRot = uti_math.vect_normalize(uti_math.vect3_prod_v(ezIn, nv))
        Mrot = uti_math.trf_rotation(vAxRot, defAng, [0]*3)[0] #Matrix of Rotation
        ex = uti_math.matr_prod(Mrot, [1,0,0])
        ey = uti_math.matr_prod(Mrot, [0,1,0])
        ez = uti_math.matr_prod(Mrot, ezIn)
        return [[tv, sv, nv], [ex, ey, ez]]

    def set_orient(self, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0): #OC21112019
        """Defines Mirror Orientation in the frame of the incident photon beam
        :param _nvx: horizontal coordinate of central normal vector
        :param _nvy: vertical coordinate of central normal vector
        :param _nvz: longitudinal coordinate of central normal vector
        :param _tvx: horizontal coordinate of central tangential vector
        :param _tvy: vertical coordinate of central tangential vector
        :param _x: horizontal position of mirror center [m]
        :param _y: vertical position of mirror center [m]
        """
        if(self.mirSub is not None): self.mirSub.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)

    def get_orient(self, _e=0): #OC18112019
        
        if(_e == 0): return self.mirSub.get_orient() #?

        nv = [self.mirSub.nvx, self.mirSub.nvy, self.mirSub.nvz]
        tv = [self.mirSub.tvx, self.mirSub.tvy, -(self.mirSub.nvx*self.mirSub.tvx + self.mirSub.nvy*self.mirSub.tvy)/self.mirSub.nvz]; 
        sv = uti_math.vect3_prod_v(nv, tv)

        #Find base vectors of the output beam (in the frame of incident beam)
        lamb = _Light_eV_mu*1.e-06/_e
        m_lamb_k0 = self.m*lamb*self.grDen*1000
        alpha = acos(-nv[2]) #?
        beta = asin(m_lamb_k0 - sin(alpha))
        beta_mi_alpha = beta - alpha
        
        ex = None; ey = None; ez = None
        angTol = 1e-07
        
        if(abs(beta_mi_alpha) <= angTol):
            ex = [-1,0,0]
            ey = [0,1,0]
            ez = [0,0,-1]
        else:
            ezIn = [0,0,1]
            vAxRot = uti_math.vect_normalize(uti_math.vect3_prod_v(ezIn, nv))
            defAng = pi + beta_mi_alpha
            Mrot = uti_math.trf_rotation(vAxRot, defAng, [0]*3)[0] #Matrix of Rotation
            ex = uti_math.matr_prod(Mrot, [1,0,0])
            ey = uti_math.matr_prod(Mrot, [0,1,0])
            ez = uti_math.matr_prod(Mrot, ezIn)
            #DEBUG
            #print('Grating Loc. Frame nv=', nv, ' ez=', ez)
            #print('norm_ex=', uti_math.vect_norm(ex), 'norm_ey=', uti_math.vect_norm(ey), 'norm_ez=', uti_math.vect_norm(ez))
            #print('Grating Mrot and [ex, ey, ez]:')
            #print(Mrot)
            #print([ex, ey, ez])
            return [[tv, sv, nv], [ex, ey, ez], Mrot] #[2] is optional transposed of [ex, ey, ez], to be used for fast space transformation calc.

        return [[tv, sv, nv], [ex, ey, ez]] #[2] is optional transposed of [ex, ey, ez], to be used for fast space transformation calc.

class SRWLOptCryst(SRWLOpt):
    """Optical Element: Ideal Crystal"""

    #def _init_(self, _d_space, _psiOr, _psiOi, _psiHr, _psiHi, _psiHBr, _psiHBi, _H1, _H2, _H3, _Tc, _Tasym, _nx, _ny, _nz, _sx, _sy, _sz, _aChi, _aPsi, _aThe):
    #def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _h1, _h2, _h3, _tc, _ang_as, _nvx, _nvy, _nvz, _tvx, _tvy):
    #def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0):
    #def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _uc=1):
    #def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx=None, _nvy=None, _nvz=None, _tvx=None, _tvy=None, _uc=1, _e_avg=0, _ang_roll=0): #OC18112019
    def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as,  _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _uc=1, _e_avg=0, _ang_roll=0): #OC18112019
        """
        :param _d_sp: (_d_space) crystal reflecting planes d-spacing (John's dA) [A]
        :param _psi0r: real part of 0-th Fourier component of crystal polarizability (John's psi0c.real) (units?)
        :param _psi0i: imaginary part of 0-th Fourier component of crystal polarizability (John's psi0c.imag) (units?)
        :param _psi_hr: (_psiHr) real part of H-th Fourier component of crystal polarizability (John's psihc.real) (units?)
        :param _psi_hi: (_psiHi) imaginary part of H-th Fourier component of crystal polarizability (John's psihc.imag) (units?)
        :param _psi_hbr: (_psiHBr:) real part of -H-th Fourier component of crystal polarizability (John's psimhc.real) (units?)
        :param _psi_hbi: (_psiHBi:) imaginary part of -H-th Fourier component of crystal polarizability (John's psimhc.imag) (units?)
        :param _tc: crystal thickness [m] (John's thicum)
        :param _ang_as: (_Tasym) asymmetry angle [rad] (John's alphdg)
        :param _nvx: horizontal coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
        :param _nvy: vertical coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
        :param _nvz: longitudinal coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
        :param _tvx: horizontal coordinate of central tangential vector (John's angles: thdg, chidg, phidg)
        :param _tvy: vertical coordinate of central tangential vector (John's angles: thdg, chidg, phidg)
        :param _uc: crystal use case: 1- Bragg Reflection, 2- Bragg Transmission (Laue cases to be added)
        :param _e_avg: average photon energy [eV] the crystal should be alligned for: if it is not 0, and _nvx, _nvy, _nvz, _tvx, _tvy are not defined, then crystal orientation will be calculated based on _uc, _e_avg and _ang_roll
        :param _ang_roll: roll angle [rad] (i.e. angle of diffraction plane angle rotation about incident beam axis) the crystal should be alligned for: it is only taken into account when _e_avg != 0
        """
        #"""
        #The Miller Inices are removed from this input (after discussion with A. Suvorov), because _d_sp already incorporates this information:
        #:param _h1: 1st index of diffraction vector (John's hMilND)
        #:param _h2: 2nd index of diffraction vector (John's kMilND)
        #:param _h3: 3rd index of diffraction vector (John's lMilND)
        #However, a member-function may be added here to calculate _d_sp from teh Miller Indices and material constant(s)
        
        #Moved to Propagation Parameters
        #:param _sx: horizontal coordinate of optical axis after crystal [m] (John's)
        #:param _sy: vertical coordinate of optical axis after crystal [m] (John's)
        #:param _sz: longitudinal coordinate of optical axis after crystal [m] (John's)

        #to go to member functions (convenience derived parameters)
        #:param _aChi: crystal roll angle (John's)
        #:param _aPsi: crystal yaw angle (John's)
        #:param _aThe: crystal theta angle (John's)
        #"""

        #DEBUG
        #print(_d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx, _nvy, _nvz, _tvx, _tvy, _uc)
        #END DEBUG
        
        self.dSp = _d_sp
        self.psi0r = _psi0r
        self.psi0i = _psi0i
        self.psiHr = _psi_hr
        self.psiHi = _psi_hi
        self.psiHbr = _psi_hbr
        self.psiHbi = _psi_hbi
        #self.h1 = _h1
        #self.h2 = _h2
        #self.h3 = _h3
        self.tc = _tc
        self.angAs = _ang_as
        self.nvx = _nvx
        self.nvy = _nvy
        self.nvz = _nvz
        self.tvx = _tvx
        self.tvy = _tvy

        self.uc = _uc #OC04092016

        self.aux_energy = _e_avg #OC17112019 #MR01082016: renamed self.energy to self.aux_energy.
        self.aux_ang_dif_pl = _ang_roll #OC17112019 #MR01082016: renamed self.ang_dif_pl to self.aux_ang_dif_pl. 
        #self.aux_energy = None  #MR01082016: renamed self.energy to self.aux_energy.
        #self.aux_ang_dif_pl = None  #MR01082016: renamed self.ang_dif_pl to self.aux_ang_dif_pl. 

        #orientIsNotDef = False if((_nvx is None) or (_nvy is None) or (_nvz is None) or (_tvx is None) or (_tvy is None)) else True #OC18112019
        #if(orientIsNotDef and (_e_avg > 0)): #OC17112019
        if(_e_avg > 0): #OC23112019
            orntData = self.find_orient(_e_avg, _ang_roll, _uc) #returns [[tv, sv, nv], [ex, ey, ez]]
            opElBaseVects = orntData[0]
            tv = opElBaseVects[0]
            nv = opElBaseVects[2]
            self.set_orient(_nvx=nv[0], _nvy=nv[1], _nvz=nv[2], _tvx=tv[0], _tvy=tv[1])
            #self.nvx = nv[0]
            #self.nvy = nv[1]
            #self.nvz = nv[2]
            #self.tvx = tv[0]
            #self.tvy = tv[1]
            #DEBUG
            #print('Beam frame vectors after Crystal:')
            #print('ex=', orntData[1][0], 'ez=', orntData[1][2])

    def set_orient(self, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0):
        """Defines Crystal Orientation in the frame of the Incident Photon beam
        :param _nvx: horizontal coordinate of normal vector
        :param _nvy: vertical coordinate of normal vector
        :param _nvz: longitudinal coordinate of normal vector
        :param _tvx: horizontal coordinate of tangential vector
        :param _tvy: vertical coordinate of tangential vector
        """
        self.nvx = _nvx
        self.nvy = _nvy
        self.nvz = _nvz
        self.tvx = _tvx
        self.tvy = _tvy

    def find_orient(self, _en, _ang_dif_pl=0, _uc=1):
        """Finds optimal crystal orientation in the input beam frame (i.e. surface normal and tangential vectors) and the orientation of the output beam frame (i.e. coordinates of the longitudinal and horizontal vectors in the input beam frame)
        :param _en: photon energy [eV]
        :param _ang_dif_pl: diffraction plane angle (0 corresponds to the vertical deflection; pi/2 to the horizontal deflection; any value in between is allowed)
        :param _uc: crystal use case: 1- Bragg Reflection, 2- Bragg Transmission (Laue cases to be added)
        :return: list of two triplets of vectors:
                out[0] is the list of 3 base vectors [tangential, saggital, normal] defining the crystal orientation
                out[1] is the list of 3 base vectors [ex, ey, ez] defining orientation of the output beam frame
                the cartesian coordinates of all these vectors are given in the frame of the input beam
        """

        #self.aux_energy = _en #MR01082016: renamed self.energy to self.aux_energy.
        #self.aux_ang_dif_pl = _ang_dif_pl #MR01082016: renamed self.ang_dif_pl to self.aux_ang_dif_pl. 
        #OC17112019 (commented-out the above)

        def prodV(_a, _b):
            return [_a[1]*_b[2] - _a[2]*_b[1], _a[2]*_b[0] - _a[0]*_b[2], _a[0]*_b[1] - _a[1]*_b[0]]
        
        def prodMV(_m, _v):
            return [_m[0][0]*_v[0] + _m[0][1]*_v[1] + _m[0][2]*_v[2],
                    _m[1][0]*_v[0] + _m[1][1]*_v[1] + _m[1][2]*_v[2],
                    _m[2][0]*_v[0] + _m[2][1]*_v[1] + _m[2][2]*_v[2]]
        
        def normV(_a):
            return sqrt(sum(n*n for n in _a)) #OC24052020
            #return sqrt(sum(n**2 for n in _a))

        #dSi = 5.43096890 # Si lattice constant (A)
        eV2wA = 12398.4193009 # energy to wavelength conversion factor 12398.41930092394
        wA = eV2wA/_en
        #Tar = math.pi*Ta/180.
        #Kh = norm(Hr)/dSi # reflection vector modulus
        kh = 1./self.dSp # because self.dSp = dSi/norm(Hr)
        hv = [0, kh*cos(self.angAs), -kh*sin(self.angAs)]
        tBr = asin(wA*kh/2)
        tKin = tBr - self.angAs # TKin = Tbr - Tar
        tKou = tBr + self.angAs # TKou = Tbr + Tar
        abs_c0 = sqrt(self.psi0r*self.psi0r + self.psi0i*self.psi0i)
        #dTref = abs(c0)*(1+math.sin(TKou)/math.sin(TKin))/2/math.sin(2*Tbr)
        dTref = 0.5*abs_c0*(1 + sin(tKou)/sin(tKin))/sin(2*tBr)
        tIn = tKin + dTref

        #Crystal orientation vectors
        nv = [0, cos(tIn), -sin(tIn)]
        tv = [0, sin(tIn), cos(tIn)]
        sv = prodV(nv, tv)

        mc = [[sv[0], nv[0], tv[0]],
              [sv[1], nv[1], tv[1]],
              [sv[2], nv[2], tv[2]]]

        #Diffracted beam frame vectors
        z1c = [sv[2], sqrt(1. - sv[2]**2 - (tv[2] + wA*hv[2])**2), tv[2] + wA*hv[2]]
        
        #DEBUG
        #print('Crystal find_orient: z1c=', z1c)
        #print('Crystal find_orient: mc=', mc)
        
        #rz = [sv[0]*z1c[0] + nv[0]*z1c[1] + tv[0]*z1c[2],
        #      sv[1]*z1c[0] + nv[1]*z1c[1] + tv[1]*z1c[2],
        #      sv[2]*z1c[0] + nv[2]*z1c[1] + tv[2]*z1c[2]]
        rz = prodMV(mc, z1c)
        
        x1c = prodV(hv, z1c)
        if sum(n**2 for n in x1c) == 0:
            x1c = prodV(nv, z1c)
        if sum(n**2 for n in x1c) == 0:
            x1c = sv
        x1c = [n/normV(x1c) for n in x1c]
        
        #rx = [sv[0]*x1c[0] + nv[0]*x1c[1] + tv[0]*x1c[2],
        #      sv[1]*x1c[0] + nv[1]*x1c[1] + tv[1]*x1c[2],
        #      sv[2]*x1c[0] + nv[2]*x1c[1] + tv[2]*x1c[2]]
        rx = prodMV(mc, x1c)
        ry = prodV(rz, rx)
        #print('ex0=',rx, 'ey0=',ry, 'ez0=',rz)

        #OC06092016
        tvNew = None; svNew = None; nvNew = None
        ex = None; ey = None; ez = None

        #The following corresponds to Bragg case (Reflection and Transmission)
        tolAng = 1.e-06
        if(abs(_ang_dif_pl) < tolAng): #case of the vertical deflection plane
            #return [[tv, sv, nv], [rx, ry, rz]]
            tvNew = tv; svNew = sv; nvNew = nv
            ex = rx; ey = ry; ez = rz
        
        else: #case of a tilted deflection plane
            cosA = cos(_ang_dif_pl)
            sinA = sin(_ang_dif_pl)
            mr = [[cosA, -sinA, 0],
                  [sinA, cosA, 0],
                  [0, 0, 1]]

            ez = prodMV(mr, rz)
            
            #Selecting "Horizontal" and "Vertical" directions of the Output beam frame
            #trying to use "minimum deviation" from the corresponding "Horizontal" and "Vertical" directions of the Input beam frame
            ezIn = [0, 0, 1]
            e1 = prodV(ez, ezIn)
            abs_e1x = abs(e1[0])
            abs_e1y = abs(e1[1])

            #ex = None; ey = None
            if(abs_e1x >= abs_e1y):
                if(e1[0] > 0): ex = e1
                else: ex = [-e1[0], -e1[1], -e1[2]]
                ex = [n/normV(ex) for n in ex]
                ey = prodV(ez, ex)
            else:
                if(e1[1] > 0): ey = e1
                else: ey = [-e1[0], -e1[1], -e1[2]]
                ey = [n/normV(ey) for n in ey]
                ex = prodV(ey, ez)
            #return [[prodMV(mr, tv), prodMV(mr, sv), prodMV(mr, nv)], [ex, ey, ez]]
            tvNew = prodMV(mr, tv) #OC06092016
            svNew = prodMV(mr, sv)
            nvNew = prodMV(mr, nv)

        #OC06092016
        #if(self.uc == 2): #Bragg Transmission
        #OC17112019
        if(_uc == 2): #Bragg Transmission
            ex = [1, 0, 0]
            ey = [0, 1, 0]
            ez = [0, 0, 1]  
        #To check this and implement other _uc cases!

        #DEBUG
        #print('Crystal find_orient:')
        #print('[tv,sv,nv]=', [tvNew, svNew, nvNew])
        return [[tvNew, svNew, nvNew], [ex, ey, ez]] #[2] can be optional transposed matrix of [ex, ey, ez], to be used for fast space transformation calc.

    def get_ang_inc(self, _e): #OC23112019
        """Estimates incidence angle for given photon energy (close to Bragg angle?)
        :param _e: photon energy [eV]
        :return: estimated incidence angle [rad]
        """        
        eV2wA = 12398.4193009 # energy to wavelength conversion factor 12398.41930092394
        wA = eV2wA/_e
        #Tar = math.pi*Ta/180.
        #Kh = norm(Hr)/dSi # reflection vector modulus
        kh = 1./self.dSp # because self.dSp = dSi/norm(Hr)
        hv = [0, kh*cos(self.angAs), -kh*sin(self.angAs)]
        tBr = asin(wA*kh/2)
        tKin = tBr - self.angAs # TKin = Tbr - Tar
        tKou = tBr + self.angAs # TKou = Tbr + Tar
        abs_c0 = sqrt(self.psi0r*self.psi0r + self.psi0i*self.psi0i)
        #dTref = abs(c0)*(1+math.sin(TKou)/math.sin(TKin))/2/math.sin(2*Tbr)
        dTref = 0.5*abs_c0*(1 + sin(tKou)/sin(tKin))/sin(2*tBr)
        return tKin + dTref

    def estim_en_fr_ang_inc(self, _ang): #OC23112019
        """Estimates photon energy from given incidence angle (in kinematic approx.)
        :param _ang: incidence angle [rad]
        :return: estimated photon energy [eV]
        """
        tBr = _ang + self.angAs
        wA = 2*self.dSp*sin(tBr)
        eV2wA = 12398.4193009 # energy to wavelength conversion factor 12398.41930092394
        return eV2wA/wA

    def get_orient(self, _e=0): #OC17112019
        """Returns data on orientation of optical element in the frame of incident beam (should be called after the optical element was completely set up)
        :return: list of two triplets of vectors:
                out[0] is the list of 3 base vectors [tangential, saggital, normal] defining the crystal orientation
                out[1] is the list of 3 base vectors [ex, ey, ez] defining orientation of the output beam frame
                the cartesian coordinates of all these vectors are given in the frame of the input beam
        """

#        def prodV(_a, _b): return [_a[1]*_b[2] - _a[2]*_b[1], _a[2]*_b[0] - _a[0]*_b[2], _a[0]*_b[1] - _a[1]*_b[0]]
#        def normV(_a): return sqrt(sum(n**2 for n in _a))
#        def prodMV(_m, _v):
#            return [_m[0][0]*_v[0] + _m[0][1]*_v[1] + _m[0][2]*_v[2],
#                    _m[1][0]*_v[0] + _m[1][1]*_v[1] + _m[1][2]*_v[2],
#                    _m[2][0]*_v[0] + _m[2][1]*_v[1] + _m[2][2]*_v[2]]

        #Crystal orientation vectors
        nv = [self.nvx, self.nvy, self.nvz]
        tv = [self.tvx, self.tvy, -(self.nvx*self.tvx + self.nvy*self.tvy)/self.nvz] #note: self.nvz < 0 for mirrors / crystals
        sv = uti_math.vect3_prod_v(nv, tv)
        exIn = [1,0,0]; eyIn = [0,1,0]; ezIn = [0,0,1]

        #print('Crystal get_orient:')
        #print('[tv,sv,nv]=', [tv, sv, nv])

        uc = 1 #Bragg Reflection use case by default
        if(hasattr(self, 'uc')):
            if(self.uc is not None): uc = self.uc

        if(uc == 2): #Bragg Transmission
            return [[tv, sv, nv], [exIn, eyIn, ezIn]] #[2] can be optional transposed matrix of [ex, ey, ez], to be used for fast space transformation calc.
        
        elif(uc == 1): #Bragg Reflection
            
            en = _e
            if(en <= 0):
                if(hasattr(self, 'aux_energy')):
                    if(self.aux_energy is not None): en = self.aux_energy

            sinAngInc = -nv[2] #?
            angInc = asin(sinAngInc)
            cosAngInc = cos(angInc)
            if(en <= 0): en = self.estim_en_fr_ang_inc(angInc)

            eV2wA = 12398.4193009 #energy to wavelength conversion factor 12398.41930092394
            wA = eV2wA/en
                
            cosAngC = cosAngInc - (wA/self.dSp)*sin(self.angAs)
            sinAngC = sqrt(1 - cosAngC*cosAngC)
            
            sinAngDef = cosAngInc*sinAngC + sinAngInc*cosAngC
            cosAngDef = cosAngInc*cosAngC - sinAngInc*sinAngC
            
            defAng = asin(sinAngDef)
            if(cosAngDef < 0): defAng = pi - defAng
        
            vAxRot = uti_math.vect_normalize(uti_math.vect3_prod_v(ezIn, nv))
            Mrot = uti_math.trf_rotation(vAxRot, defAng, [0]*3)[0] #Matrix of Rotation
            ex = uti_math.matr_prod(Mrot, exIn)
            ey = uti_math.matr_prod(Mrot, eyIn)
            ez = uti_math.matr_prod(Mrot, ezIn)
            
            return [[tv, sv, nv], [ex, ey, ez], Mrot] #[2] is optional transposed of [ex, ey, ez], to be used for fast space transformation calc.
          
        #elif(uc == 3): #Laue
        #elif(uc == 4): #Laue Transmission

        return [[tv, sv, nv], [exIn, eyIn, ezIn]]

class SRWLOptC(SRWLOpt):
    """Optical Element: Container"""
    
    def __init__(self, _arOpt=None, _arProp=None):
        """
        :param _arOpt: optical element structures list (or array)
        :param _arProp: list of lists of propagation parameters to be used for each individual optical element
            Each element _arProp[i] is a list in which elements mean:
            [0]: Auto-Resize (1) or not (0) Before propagation
            [1]: Auto-Resize (1) or not (0) After propagation
            [2]: Relative Precision for propagation with Auto-Resizing (1. is nominal)
            [3]: Allow (1) or not (0) for semi-analytical treatment of the quadratic (leading) phase terms at the propagation
            [4]: Do any Resizing on Fourier side, using FFT, (1) or not (0)
            [5]: Horizontal Range modification factor at Resizing (1. means no modification)
            [6]: Horizontal Resolution modification factor at Resizing (1. means no modification)
            [7]: Vertical Range modification factor at Resizing (1. means no modification)
            [8]: Vertical Resolution modification factor at Resizing (1. means no modification)
            [9]: Optional: Type of wavefront Shift before Resizing (vs which coordinates; to be implemented)
            [10]: Optional: New Horizontal wavefront Center position after Shift (to be implemented)
            [11]: Optional: New Vertical wavefront Center position after Shift (to be implemented)
            [12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
            [13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
            [14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
            [15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
            [16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate
        """
        self.arOpt = _arOpt #optical element structures array
        if(_arOpt is None):
            self.arOpt = []
        self.arProp = _arProp #list of lists of propagation parameters to be used for individual optical elements
        if(_arProp is None):
            self.arProp = []
            
    def allocate(self, _nElem):
        self.arOpt = [SRWLOpt()]*_nElem
        self.arProp = [[0]*17]*_nElem

    def append_drift(self, _len): #OC28042018
        """Appends drift space to the end of the Container
        :param _len: length [m]
        """
        if(_len == 0.): return

        nOE = len(self.arOpt)
        lastOptEl = self.arOpt[nOE - 1]

        if(isinstance(lastOptEl, SRWLOptD)): #just change the length of the last drift
            lastOptEl.L += _len
            #print('Drift length updated to:', lastOptEl.L)
        else:
            self.arOpt.append(SRWLOptD(_len))
            ppDr = [0, 0, 1.0, 1, 0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #to improve
            nPP = len(self.arProp)
            nOE += 1
            if(nOE == nPP):
                ppPost = copy(self.arProp[nPP - 1])
                self.arProp[nPP - 1] = ppDr
                self.arProp.append(ppPost)
            elif(nOE > nPP):
                self.arProp.append(ppDr)

    def get_orient_lab_fr(self, _e=0, _r=0, _v_op_ax=[0,0,1]):
        """Returns Cartesian coordinates of all optical elements' center positions and base vectors (t, s, n)
        :param _e: photon energy [eV] optical scheme shoule be represented for (required for crystals and gratings)
        :param _r: distance to first optical element (along beam axis) [m] or list of Cartesian coordinates of the center of the first optical element
        :param _v_op_ax: Cartesian coordinates of the initial beam axis vector in the lab frame [m]
        """
        resData = []
        nElem = len(self.arOpt)
        if(nElem <= 0): return resData

        opAxVect = copy(_v_op_ax) if(_v_op_ax is not None) else [0,0,1] #Beam axis

        P = _r if(isinstance(_r, list) or isinstance(_r, array)) else uti_math.vect_mult(copy(opAxVect), _r) #Optical elememt position
        #P = uti_math.vect_mult(copy(opAxVect), _r) #Optical elememt position
        
        bmFrTrfMatr = [[1,0,0],[0,1,0],[0,0,1]] #Beam frame transformation matrix
        lenDrift = 0
        
        for i in range(nElem):
            opEl = self.arOpt[i]

            if(isinstance(opEl, SRWLOptD)):
                lenDrift += opEl.L
                for k in range(3): P[k] += lenDrift*opAxVect[k] #Updating (translating along beam axis) position of next optical element 
            
            elif(isinstance(opEl, SRWLOptC)): 
                resData.append(opEl.get_orient_lab_fr(_e=_e, _r=P, _v_op_ax=opAxVect)) #?
                #resData.append(opEl.get_orient_lab_fr(_e=_e, _r=lenDrift, _v_op_ax=opAxVect)) #Perhaps initial point should be submitted instead of lenDrift?
                lenDrift = 0

            else:
                #Get tx,ty,tz,sx,sy,sz,nx,ny,nz from opEl
                opElOrntData = opEl.get_orient(_e)
                [tv, sv, nv] = opElOrntData[0]
                
                tv = uti_math.matr_prod(bmFrTrfMatr, tv) #Transforming coordinates of tv, sv, nv vectors to Lab frame 
                sv = uti_math.matr_prod(bmFrTrfMatr, sv)
                nv = uti_math.matr_prod(bmFrTrfMatr, nv)

                opElClassName = opEl.__class__.__name__ #Does it work in Py 2.7?
                opElID = opElClassName[7:] #To get rid of "SRWLOpt"
                curOpElData = [opElID, P[0],P[1],P[2], tv[0],tv[1],tv[2], sv[0],sv[1],sv[2], nv[0],nv[1],nv[2]] #to become [x,y,z,tx,ty,tz,sx,sy,sz,nx,ny,nz]
                resData.append(curOpElData)
                
                #opElExEyEzMatr = opElOrntData[1]
                opElTrfMatr = uti_math.matr_transp(opElOrntData[1]) if(len(opElOrntData) < 3) else opElOrntData[2]
                opAxVect = uti_math.matr_prod(opElTrfMatr, opAxVect) #Updating beam axis vector

                #DEBUG
                #print('Transformation Matrix Before Update:')
                #print(bmFrTrfMatr)
                #print('Transformation Matrix of Individual Opt. Elem.:')
                #print(opElOrntData[1])
                
                bmFrTrfMatr = uti_math.matr_prod(opElTrfMatr, bmFrTrfMatr) #Updating transformation matrix by multiplying it by new opt. elem. matrix from left
                lenDrift = 0
                
                #DEBUG
                #print('Transformation Matrix After Update:')
                #print(bmFrTrfMatr)
                #print('Beam Frame Base Vector Coords. after ', opElID, ':')
                #print('ex=', bmFrTrfMatr[0], '  ey=', bmFrTrfMatr[1], '  ez=', bmFrTrfMatr[2])
                #print('norm_ex=', uti_math.vect_norm(bmFrTrfMatr[0]), 'norm_ey=', uti_math.vect_norm(bmFrTrfMatr[1]), 'norm_ez=', uti_math.vect_norm(bmFrTrfMatr[2]))

        return resData

    def randomize(self):
        """Overrides SRWLOpt.randomize(); randomizes parameters of optical elements in the container
        """
        if(not hasattr(self, 'arOpt')): return
        if(self.arOpt is None): return
        if(not isinstance(self.arOpt, list)): return
        lenArOpt = len(self.arOpt)
        if(lenArOpt <= 0): return
        for i in range(lenArOpt): self.arOpt[i].randomize()

#****************************************************************************
#****************************************************************************
#Setup some transmission-type optical elements
#****************************************************************************
#****************************************************************************
#def srwl_opt_setup_CRL(_foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n, _wall_thick, _xc, _yc, _void_cen_rad=None, _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
def srwl_opt_setup_CRL(_foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n, _wall_thick, _xc, _yc, _void_cen_rad=None, _e_start=0, _e_fin=0, _nx=1001, _ny=1001, _ang_rot_ex=0, _ang_rot_ey=0): #, _ang_rot_ez=0): #OC15042018
    """
    Setup Transmission type Optical Element which simulates Compound Refractive Lens (CRL)
    :param _foc_plane: plane of focusing: 1- horizontal, 2- vertical, 3- both
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _shape: 1- parabolic, 2- circular (spherical)
    :param _apert_h: horizontal aperture size [m]
    :param _apert_v: vertical aperture size [m]
    :param _r_min: radius (on tip of parabola for parabolic shape) [m]
    :param _n: number of lenses (/"holes")
    :param _wall_thick: min. wall thickness between "holes" [m]
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :param _void_cen_rad: flat array/list of void center coordinates and radii: [x1, y1, r1, x2, y2, r2,...] 
    :param _e_start: initial photon energy
    :param _e_fin: final photon energy
    :param _nx: number of points vs horizontal position to represent the transmission element
    :param _ny: number of points vs vertical position to represent the transmission element
    :param _ang_rot_ex: angle [rad] of CRL rotation about horizontal axis
    :param _ang_rot_ey: angle [rad] of CRL rotation about vertical axis
    :param _ang_rot_ez: angle [rad] of CRL rotation about longitudinal axis
    :return: transmission (SRWLOptT) type optical element which simulates CRL
    """

    input_parms = { #MR26022016: Added all input parameters to include in return object:
        "type": "crl",
        "focalPlane": _foc_plane,
        "refractiveIndex": _delta,
        "attenuationLength": _atten_len,
        "shape": _shape,
        "horizontalApertureSize": _apert_h,
        "verticalApertureSize": _apert_v,
        "radius": _r_min,
        "numberOfLenses": _n,
        "wallThickness": _wall_thick,
        "horizontalCenterCoordinate": _xc,
        "verticalCenterCoordinate": _yc,
        "voidCenterCoordinates": _void_cen_rad,
        "initialPhotonEnergy": _e_start,
        "finalPhotonPnergy": _e_fin,
        "horizontalPoints": _nx,
        "verticalPoints": _ny,
    }

    def ray_path_in_one_CRL(_x, _y, _foc_plane, _shape, _half_apert, _r_min, _wall_thick): #CRL is always centered
        rE2 = 0
        if((_foc_plane == 1) or (_foc_plane == 3)): #focusing in horizontal plane
            rE2 += _x*_x
        if((_foc_plane == 2) or (_foc_plane == 3)): #focusing in vertical or in both planes
            rE2 += _y*_y
            
        halfApE2 = _half_apert*_half_apert
        
        sectLen = 0
        if(_shape == 1): #parabolic
            a = 1./_r_min
            sectLen = _wall_thick + a*halfApE2
            if(rE2 < halfApE2):
                return _wall_thick + a*rE2
        elif(_shape == 2): #circular (or spherical)
            radE2 = _r_min*_r_min
            sectLen = _wall_thick + 2*_r_min
            if(_half_apert < _r_min):
                sectLen = _wall_thick + 2*(_r_min - sqrt(radE2 - halfApE2))
                if(rE2 < halfApE2):
                    #return sectLen - 2*sqrt(radE2 - rE2)
                    return _wall_thick + 2*(_r_min - sqrt(radE2 - rE2)) #OC22042018
            elif(rE2 < radE2):
                return sectLen - 2*sqrt(radE2 - rE2)

        return sectLen

    #OC20042018
    def ray_path_in_one_CRL_rot(_x, _y, _foc_plane, _shape, _half_apert, _r_min, _wall_thick, _kx, _ky): #CRL is centered, but rotated about vert. / hor  axes
        rE2 = 0
        if((_foc_plane == 1) or (_foc_plane == 3)): #focusing in horizontal plane
            rE2 += _x*_x
        if((_foc_plane == 2) or (_foc_plane == 3)): #focusing in vertical or in both planes
            rE2 += _y*_y
            
        halfApE2 = _half_apert*_half_apert

        kx2 = _kx*_kx
        ky2 = _ky*_ky
        tiltFact = sqrt(1. + kx2 + ky2)
        sectLen = 0
        if(_shape == 1): #parabolic
            a = 1./_r_min
            sectLen = (_wall_thick + a*halfApE2)*tiltFact
            if(rE2 < halfApE2):
                if(_foc_plane == 1):
                    if(_kx == 0): return (_wall_thick + a*rE2)*tiltFact
                    else:
                        a_kx = a*_kx
                        a_kx2 = a_kx*_kx
                        a_kx2_w = a_kx2*_wall_thick
                        two_a_kx_x = 2*a_kx*_x
                        return tiltFact*(2 - sqrt(1 - a_kx2_w + two_a_kx_x) - sqrt(1 - a_kx2_w - two_a_kx_x))/a_kx2
                elif(_foc_plane == 2):
                    if(_ky == 0): return (_wall_thick + a*rE2)*tiltFact
                    else:
                        a_ky = a*_ky
                        a_ky2 = a_ky*_ky
                        a_ky2_w = a_ky2*_wall_thick
                        two_a_ky_y = 2*a_ky*_y
                        return tiltFact*(2 - sqrt(1 - a_ky2_w + two_a_ky_y) - sqrt(1 - a_ky2_w - two_a_ky_y))/a_ky2
                else:
                    a_k2 = a*(kx2 + ky2)
                    buf0 = _ky*_x - _kx*_y
                    buf1 = 1 - a*a*buf0*buf0 - w*a_k2
                    buf2 = 2*a*(_kx*_x + _ky*_y)
                    return tiltFact*(2 - sqrt(buf1 - buf2) - sqrt(buf1 + buf2))/a_k2

        elif(_shape == 2): #circular (or spherical)

            #This takes into account circular CRL tilt / rotation only approximately!
            radE2 = _r_min*_r_min
            sectLen = (_wall_thick + 2*_r_min)*tiltFact
            
            if(_half_apert < _r_min):
                sectLen = _wall_thick + 2*(_r_min - sqrt(radE2 - halfApE2))
                if(rE2 < halfApE2):
                    #return sectLen - 2*sqrt(radE2 - rE2)
                    return (_wall_thick + 2*(_r_min - sqrt(radE2 - rE2)))*tiltFact #OC22042018
                    
            elif(rE2 < radE2):
                return (sectLen - 2*sqrt(radE2 - rE2))*tiltFact

        #return sectLen
        return sectLen*tiltFact

    def ray_path_in_spheres(_x, _y, _void_cen_rad):
        n = int(round(len(_void_cen_rad)/3))
        sumPath = 0.
        for i in range(n):
            i3 = i*3
            dx = _x - _void_cen_rad[i3]
            dy = _y - _void_cen_rad[i3 + 1]
            rVoid = _void_cen_rad[i3 + 2]
            uE2 = dx*dx + dy*dy
            rVoidE2 = rVoid*rVoid
            if(uE2 < rVoidE2):
                sumPath += 2*sqrt(rVoidE2 - uE2)
                #print('Void crossed:', dx, dy, rVoid, sumPath)
        return sumPath

    #foc_len = (0.5*_r_min/(_n*_delta))
    #print('Optical Element Setup: CRL Focal Length:', foc_len, 'm')

    #fx = 1e+23
    #fy = 1e+23
    #if(_foc_plane != 1):
    #    fy = foc_len
    #if(_foc_plane != 2):
    #    fx = foc_len

    rx = _apert_h*1.1
    ry = _apert_v*1.1
    nx = _nx #1001
    ny = _ny #1001
    ne = 1
    arDelta = [0]
    arAttenLen = [0]

    #if(((type(_delta).__name__ == 'list') or (type(_delta).__name__ == 'array')) and ((type(_atten_len).__name__ == 'list') or (type(_atten_len).__name__ == 'array'))):
    if(isinstance(_delta, list) or isinstance(_delta, array)) and (isinstance(_atten_len, list) or isinstance(_atten_len, array)):
        ne = len(_delta)
        ne1 = len(_atten_len)
        if(ne > ne1):
            ne = ne1
        arDelta = _delta
        arAttenLen = _atten_len
    else:
        arDelta[0] = _delta
        arAttenLen[0] = _atten_len

    trM = None #OC15042018
    vZero = [0,0,0]
    if(_ang_rot_ex != 0):
        trRot = uti_math.trf_rotation([1,0,0], _ang_rot_ex, vZero)
        trM = trRot[0]
    if(_ang_rot_ey != 0):
        trRot = uti_math.trf_rotation([0,1,0], _ang_rot_ey, vZero)
        if(trM is not None): trM = uti_math.matr_prod(trRot[0], trM)
        else: trM = trRot[0]
        
    kxRay = 0; kyRay = 0
    #cosecAngX = 1; cosecAngY = 1
    #fact_f = 1
    fact_opt_path = 1
    kTol = 1.e-07 #to tune
    if(trM is not None):
        trMi = uti_math.matr_3x3_inv(trM)
        vRay = uti_math.matr_prod(trMi, [0,0,1])
        kxRay = vRay[0]/vRay[2]
        if(abs(kxRay) < kTol): kxRay = 0
        kyRay = vRay[1]/vRay[2]
        if(abs(kyRay) < kTol): kyRay = 0

        fact_opt_path = sqrt(1. + kxRay*kxRay + kyRay*kyRay) #optical path factor (due to CRL rotations)

        #if(kxRay != 0): cosecAngX = 1/vRay[0]
        #if(kyRay != 0): cosecAngY = 1/vRay[1]
        #fact_f = 1./sqrt(1. + kxRay*kxRay + kyRay*kyRay) #focal length factor (due to CRL rotations)

    #foc_len = 0.5*_r_min/(_n*arDelta[int(0.5*ne)])
    foc_len = 0.5*_r_min/(_n*arDelta[int(0.5*ne)]*fact_opt_path) #OC19042018
    if srwl_uti_proc_is_master(): print('Optical Element Setup: CRL Focal Length:', foc_len, 'm') #OC27102021
    #print('Optical Element Setup: CRL Focal Length:', foc_len, 'm')

    fx = 1e+23
    fy = 1e+23
    if(_foc_plane != 1):
        fy = foc_len
    if(_foc_plane != 2):
        fx = foc_len
    
    opT = SRWLOptT(nx, ny, rx, ry, None, 1, fx, fy, _xc, _yc, ne, _e_start, _e_fin)

    halfApert = 0.5*_apert_h
    halfApertV = 0.5*_apert_v

    if(_foc_plane == 2): #1D lens, vertical is focusing plane
        halfApert = halfApertV
    elif(_foc_plane == 3): #2D lens
        if(halfApert > halfApertV):
            halfApert = halfApertV
    
    hx = rx/(nx - 1)
    hy = ry/(ny - 1)

    fact_opt_path_n = _n #OC19042018
    complexRotCase = False
    if((kxRay != 0) or (kyRay != 0)):
        if(((kxRay != 0) and (_foc_plane == 2)) or ((kyRay != 0) and (_foc_plane == 1))): fact_opt_path_n = fact_opt_path*_n
        else: complexRotCase = True
    
    #Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    ofst = 0
    y = -0.5*ry #CRL is always centered on the grid, however grid can be shifted
    for iy in range(ny):
        x = -0.5*rx
        for ix in range(nx):
            #pathInBody = _n*ray_path_in_one_CRL(x, y, _foc_plane, _shape, halfApert, _r_min, _wall_thick)
            #OC20042018
            if complexRotCase:
                pathInBody = fact_opt_path_n*ray_path_in_one_CRL_rot(x, y, _foc_plane, _shape, halfApert, _r_min, _wall_thick, kxRay, kyRay)
            else:
                pathInBody = fact_opt_path_n*ray_path_in_one_CRL(x, y, _foc_plane, _shape, halfApert, _r_min, _wall_thick)
    
            if(_void_cen_rad is not None): #eventually subtract path in voids
                pathInBody -= ray_path_in_spheres(x, y, _void_cen_rad)

            for ie in range(ne):
                opT.arTr[ofst] = exp(-0.5*pathInBody/arAttenLen[ie]) #amplitude transmission
                opT.arTr[ofst + 1] = -arDelta[ie]*pathInBody #optical path difference
                ofst += 2
            x += hx
        y += hy

    opT.input_parms = input_parms #MR16022016
    return opT

#****************************************************************************
def srwl_opt_setup_saw_tooth_lens(_delta, _atten_len, _dxdy, _thick, _ang_wedge=0, _ang_rot=0, _per_x=0, _per_y=0, _hole_nx=1, _hole_ny=1, _xc=0, _yc=0, _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
    """
    Setup Transmission type Optical Element which simulates trapethoidal shape etched mask (requested by K. Kaznatcheev)
    :param _delta: refractive index decrement
    :param _atten_len: attenuation length [m]
    :param _dxdy: array of transverse dimensions (pairs) of base holes [m]
    :param _ang_wedge: tooth profile wedge angle [rad]
    :param _ang_rot: rotation angle [rad]
    :param _per_x: period of hole / lens pattern in the horizontal direction [m]
    :param _per_y: period of hole / lens pattern in the vertical direction [m]
    :param _hole_nx: number of holes / lenses in the horizontal direction [m]
    :param _hole_ny: number of holes / lenses in the vertical direction [m]
    :param _xc: horizontal coordinate of pattern center [m]
    :param _yc: vertical coordinate of pattern center [m]
    :param _e_start: initial photon energy [eV]
    :param _e_fin: final photon energy [eV]
    :param _nx: number of points vs horizontal position to represent the transmission element
    :param _ny: number of points vs vertical position to represent the transmission element
    :return: transmission (SRWLOptT) type optical element which simulates Cylindrical Fiber
    """

    #input_parms = {
    #    "type": "crl",
    #    "focalPlane": _foc_plane,
    #    "refractiveIndex": _delta,
    #    "attenuationLength": _atten_len,
    #    "shape": _shape,
    #    "horizontalApertureSize": _apert_h,
    #    "verticalApertureSize": _apert_v,
    #    "radius": _r_min,
    #    "numberOfLenses": _n,
    #    "wallThickness": _wall_thick,
    #    "horizontalCenterCoordinate": _xc,
    #    "verticalCenterCoordinate": _yc,
    #    "voidCenterCoordinates": _void_cen_rad,
    #    "initialPhotonEnergy": _e_start,
    #    "finalPhotonPnergy": _e_fin,
    #    "horizontalPoints": _nx,
    #    "verticalPoints": _ny,
    #}

    ar_dxdy = None
    if(isinstance(_dxdy, array) or isinstance(_dxdy, list)):
        len_dxdy = len(_dxdy)
        dxdy0 = _dxdy[0]
        if(isinstance(dxdy0, array) or isinstance(dxdy0, list)):
            len_dxdy0 = len(dxdy0)
            if(len_dxdy0 > 2): raise Exception("Incorrect definition of hole dimensions in trapethoidal shape etched mask")
            elif(len_dxdy0 == 2): ar_dxdy = _dxdy
            elif(len_dxdy0 == 1):
                ar_dxdy = [[0,0]]*len_dxdy
                for i in range(len_dxdy):
                    d = _dxdy[i][0]
                    ar_dxdy[i][0] = d
                    ar_dxdy[i][1] = d
        else:
            if(len_dxdy == 2): ar_dxdy = [[_dxdy[0], _dxdy[1]]]
            else:
                ar_dxdy = [[0,0]]*len_dxdy
                for i in range(len_dxdy):
                    d = _dxdy[i]
                    ar_dxdy[i] = [d, d]
                    #ar_dxdy[i][1] = d
                    #TEST
                    #print(d)
                    #print(ar_dxdy[i])
    else:
        ar_dxdy = [[_dxdy, _dxdy]]*len_dxdy

    #TEST
    #print(ar_dxdy)
    
    tanAng = tan(_ang_wedge)
    rExt1 = _thick/tanAng
    rExt = 2.*rExt1

    nLayers = len(ar_dxdy)
    ar_half_rx1 = [0]*nLayers
    ar_half_ry1 = [0]*nLayers
    
    rx = 0; ry = 0
    for i in range(nLayers):
        cur_dxdy = ar_dxdy[i]
        cur_rx = cur_dxdy[0]
        if(rx < cur_rx): rx = cur_rx
        cur_ry = cur_dxdy[1]
        if(ry < cur_ry): ry = cur_ry
        ar_half_rx1[i] = 0.5*(cur_rx + rExt)
        ar_half_ry1[i] = 0.5*(cur_ry + rExt)
  
    rxTot = rx + rExt*1.2;
    ryTot = ry + rExt*1.2;

    if(_hole_nx > 1): rxTot += (_hole_nx - 1)*_per_x
    if(_hole_ny > 1): ryTot += (_hole_ny - 1)*_per_y

    rx1 = rx + rExt
    ry1 = ry + rExt
    half_rx1 = 0.5*rx1
    half_ry1 = 0.5*ry1

    #TEST
    #print(rx, rExt, rExt1)

    ne = 1
    arDelta = [0]
    arAttenLen = [0]
    if(isinstance(_delta, list) or isinstance(_delta, array)) and (isinstance(_atten_len, list) or isinstance(_atten_len, array)):
        ne = len(_delta)
        ne1 = len(_atten_len)
        if(ne > ne1): ne = ne1
        arDelta = _delta
        arAttenLen = _atten_len
    else:
        arDelta[0] = _delta
        arAttenLen[0] = _atten_len

    baseAmpTr = exp(-0.5*_thick/arAttenLen[0]) #amplitude transmission
    baseOptPathDif = -arDelta[0]*_thick #opt. path dif.
    opT = SRWLOptT(_nx, _ny, rxTot, ryTot, _arTr=None, _extTr=1, _Fx=1e+23, _Fy=1e+23, _x=_xc, _y=_yc, _ne=ne, _eStart=_e_start, _eFin=_e_fin, _alloc_base=[baseAmpTr, baseOptPathDif])
    arTr = opT.arTr

    nxny = _nx*_ny
    two_ne = 2*ne
    for ie in range(1, ne):
        baseAmpTr = exp(-0.5*_thick/arAttenLen[ie]) #amplitude transmission
        baseOptPathDif = -arDelta[ie]*_thick #opt. path dif.
        ofst = ie*2
        for ixy in range(nxny):
            arTr[ofst] = baseAmpTr
            arTr[ofst + 1] = baseOptPathDif
            ofst += two_ne

    xStep = rxTot/(_nx - 1)
    xStart = _xc - 0.5*rxTot
    yStep = ryTot/(_ny - 1)
    yStart = _yc - 0.5*ryTot
    
    xcMin = _xc - 0.5*(_hole_nx - 1)*_per_x
    ycMin = _yc - 0.5*(_hole_ny - 1)*_per_y

    perIX = 2*ne
    perIY = perIX*_nx

    #TEST
    #print(ar_half_rx1)
    #print(ar_half_ry1)
    
    ycCur = ycMin
    for ihy in range(_hole_ny):
        xcCur = xcMin
        for ihx in range(_hole_nx):

            yStartCur = ycCur - half_ry1
            yEndCur = ycCur + half_ry1

            diyStartCur = (yStartCur - yStart)/yStep
            iyStartCur = int(diyStartCur)
            if((diyStartCur - iyStartCur) > 1.e-06): iyStartCur += 1
            #TEST
            #print(diyStartCur, iyStartCur, diyStartCur - iyStartCur)
            
            if(iyStartCur < 0): iyStartCur = 0
            elif(iyStartCur >= _ny): iyStartCur = _ny - 1

            iyEndCur = int((yEndCur - yStart)/yStep) + 1
            if(iyEndCur < 0): iyEndCur = 0
            elif(iyEndCur > _ny): iyStartCur = _ny

            xStartCur = xcCur - half_rx1
            xEndCur = xcCur + half_rx1

            dixStartCur = (xStartCur - xStart)/xStep
            ixStartCur = int(dixStartCur)
            if((dixStartCur - ixStartCur) > 1.e-06): ixStartCur += 1
            #TEST
            #print(dixStartCur, ixStartCur, dixStartCur - ixStartCur)
            
            if(ixStartCur < 0): ixStartCur = 0
            elif(ixStartCur >= _nx): ixStartCur = _nx - 1
            ixEndCur = int((xEndCur - xStart)/xStep) + 1
            if(ixEndCur < 0): ixEndCur = 0
            elif(ixEndCur > _nx): ixStartCur = _nx

            #TEST
            #print(iyStartCur, iyEndCur)
            #print(ixStartCur, ixEndCur)

            for iy in range(iyStartCur, iyEndCur):
                y = yStart + yStep*iy

                for ix in range(ixStartCur, ixEndCur):
                    x = xStart + xStep*ix

                    pathInBody = 0
                    for i in range(nLayers):
                        cur_half_ry1 = ar_half_ry1[i]
                        yStartCurLayer = ycCur - cur_half_ry1
                        yEndCurLayer = ycCur + cur_half_ry1
                        yStartCurLayer1 = yStartCurLayer + rExt1
                        yEndCurLayer1 = yEndCurLayer - rExt1

                        cur_half_rx1 = ar_half_rx1[i]
                        xStartCurLayer = xcCur - cur_half_rx1
                        xEndCurLayer = xcCur + cur_half_rx1
                        xStartCurLayer1 = xStartCurLayer + rExt1
                        xEndCurLayer1 = xEndCurLayer - rExt1

                        #TEST
                        #print(yStartCurLayer, yStartCurLayer1, xStartCurLayer1, xEndCurLayer1)
                        #if((ix == 72) and (iy == 72)):
                        #    print(yStartCurLayer, y, yStartCurLayer1)
                        #    print(xStartCurLayer, x, xStartCurLayer1)

                        if(yStartCurLayer <= y <= yStartCurLayer1):
                            if(xStartCurLayer <= x <= xStartCurLayer1):
                                b = yStartCurLayer - xStartCurLayer
                                if(y > x + b): pathInBody += (xStartCurLayer1 - x)*tanAng
                                else: pathInBody += (yStartCurLayer1 - y)*tanAng
                            elif(xStartCurLayer1 < x < xEndCurLayer1):
                                pathInBody += (yStartCurLayer1 - y)*tanAng
                                #TEST
                                #print(y, pathInBody)

                            elif(xEndCurLayer1 <= x <= xEndCurLayer):
                                b = yStartCurLayer + xEndCurLayer
                                if(y < b - x): pathInBody += (yStartCurLayer1 - y)*tanAng
                                else: pathInBody += (x - xEndCurLayer1)*tanAng
                                
                        elif(yStartCurLayer1 < y < yEndCurLayer1):
                            if(xStartCurLayer <= x <= xStartCurLayer1): pathInBody += (xStartCurLayer1 - x)*tanAng
                            elif(xEndCurLayer1 <= x <= xEndCurLayer): pathInBody += (x - xEndCurLayer1)*tanAng

                        elif(yEndCurLayer1 <= y <= yEndCurLayer):
                            if(xStartCurLayer <= x <= xStartCurLayer1):
                                b = yEndCurLayer + xStartCurLayer
                                if(y < b - x): pathInBody += (xStartCurLayer1 - x)*tanAng
                                else: pathInBody += (y - yEndCurLayer1)*tanAng
                            elif(xStartCurLayer1 < x < xEndCurLayer1): pathInBody += (y - yEndCurLayer1)*tanAng
                            elif(xEndCurLayer1 <= x <= xEndCurLayer):
                                b = yEndCurLayer - xEndCurLayer
                                if(y > x + b): pathInBody += (y - yEndCurLayer1)*tanAng
                                else: pathInBody += (x - xEndCurLayer1)*tanAng

                    #if(pathInBody > 0):
                    ofst = iy*perIY + ix*perIX
                    for ie in range(ne):
                        opT.arTr[ofst] = exp(-0.5*pathInBody/arAttenLen[ie]) #amplitude transmission
                        opT.arTr[ofst + 1] = -arDelta[ie]*pathInBody #optical path difference
                        ofst += 2
                
            xcCur += _per_x
        ycCur += _per_y

    #opT.input_parms = input_parms
    return opT

#****************************************************************************
def srwl_opt_setup_cyl_fiber(_foc_plane, _delta_ext, _delta_core, _atten_len_ext, _atten_len_core, _diam_ext, _diam_core, _xc, _yc):
    """
    Setup Transmission type Optical Element which simulates Cylindrical Fiber
    :param _foc_plane: plane of focusing: 1- horizontal (i.e. fiber is parallel to vertical axis), 2- vertical (i.e. fiber is parallel to horizontal axis)
    :param _delta_ext: refractive index decrement of external layer
    :param _delta_core: refractive index decrement of core
    :param _atten_len_ext: attenuation length [m] of external layer
    :param _atten_len_core: attenuation length [m] of core
    :param _diam_ext: diameter [m] of external layer
    :param _diam_core: diameter [m] of core
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :return: transmission (SRWLOptT) type optical element which simulates Cylindrical Fiber
    """

    input_parms = { #MR26022016: Added all input parameters to include in return object:
        "type": "cyl_fiber",
        "focalPlane": _foc_plane,
        "externalRefractiveIndex": _delta_ext,
        "coreRefractiveIndex": _delta_core,
        "externalAttenuationLength": _atten_len_ext,
        "coreAttenuationLength": _atten_len_core,
        "externalDiameter": _diam_ext,
        "coreDiameter": _diam_core,
        "horizontalCenterPosition": _xc,
        "verticalCenterPosition": _yc,
    }

    def ray_path_in_cyl(_dx, _diam):
        r = 0.5*_diam
        pathInCyl = 0
        if((_dx > -r) and (_dx < r)):
            pathInCyl = 2*sqrt(r*r - _dx*_dx)
        return pathInCyl

    ne = 1
    nx = 101
    ny = 1001
    rx = 10e-03
    ry = _diam_ext*1.2
    if(_foc_plane == 1): #focusing plane is horizontal
        nx = 1001
        ny = 101
        rx = _diam_ext*1.2
        ry = 10e-03

    opT = SRWLOptT(nx, ny, rx, ry, None, 1, 1e+23, 1e+23, _xc, _yc)

    hx = rx/(nx - 1)
    hy = ry/(ny - 1)
    ofst = 0
    pathInExt = 0
    pathInCore = 0

    if(_foc_plane == 2): #focusing plane is vertical
        y = -0.5*ry #cylinder is always centered on the grid, however grid can be shifted
        for iy in range(ny):
            pathInExt = 0; pathInCore = 0
            if(_diam_core > 0):
                pathInCore = ray_path_in_cyl(y, _diam_core)
            pathInExt = ray_path_in_cyl(y, _diam_ext) - pathInCore
            argAtten = -0.5*pathInExt/_atten_len_ext
            if(_atten_len_core > 0):
                argAtten -= 0.5*pathInCore/_atten_len_core
            ampTr = exp(argAtten) #amplitude transmission
            optPathDif = -_delta_ext*pathInExt - _delta_core*pathInCore #optical path difference
            for ix in range(nx):                    
                opT.arTr[ofst] = ampTr #amplitude transmission
                opT.arTr[ofst + 1] = optPathDif #optical path difference
                ofst += 2
            y += hy
    else: #focusing plane is horizontal
        perY = 2*nx
        x = -0.5*rx #cylinder is always centered on the grid, however grid can be shifted
        for ix in range(nx):
            pathInExt = 0; pathInCore = 0
            if(_diam_core > 0):
                pathInCore = ray_path_in_cyl(x, _diam_core)
            pathInExt = ray_path_in_cyl(x, _diam_ext) - pathInCore
            argAtten = -0.5*pathInExt/_atten_len_ext
            if(_atten_len_core > 0):
                argAtten -= 0.5*pathInCore/_atten_len_core
            ampTr = exp(argAtten) #amplitude transmission
            optPathDif = -_delta_ext*pathInExt - _delta_core*pathInCore #optical path difference
            ix2 = ix*2
            for iy in range(ny):
                ofst = iy*perY + ix2
                opT.arTr[ofst] = ampTr #amplitude transmission
                opT.arTr[ofst + 1] = optPathDif #optical path difference
            x += hx

    opT.input_parms = input_parms  #MR26022016
    return opT

#****************************************************************************
#OC: To rename "mask" to something more meaningful, using general physical/technical terms
#OC27012019: This probably needs to be fixed / re-programmed
#OC: Failed attempt to use it:
#opG = srwl_opt_setup_mask(_delta=1e-06, _atten_len=30e-06, _thick=1.e-03,
#                          _hx=1.e-06, _hy=1.e-06, _pitch_x=1000e-06, _pitch_y=20.e-06, _mask_Nx=1000, _mask_Ny=1000,
#                          _grid_nx=1, _grid_ny=100, _grid_sh=1, _grid_dx=1.e-03, _grid_dy=8.e-06, _grid_angle=0, _mask_x0=0, _mask_y0=0)
def srwl_opt_setup_mask(_delta, _atten_len, _thick,
                        _hx, _hy, _pitch_x, _pitch_y, _mask_Nx, _mask_Ny,
                        _grid_nx, _grid_ny, _grid_sh, _grid_dx, _grid_dy=0, _grid_angle=0, _mask_x0=0, _mask_y0=0):
    """Setup Transmission type Optical Element which simulates a mask array for at-wavelength metrology.

    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _thick: thickness of mask [m]
    :param _hx: sampling interval in x-direction [m]
    :param _hy: sampling interval in y-direction [m]
    :param _pitch_x: grid pitch in x-direction [m]
    :param _pitch_y: grid pitch in y-direction [m]
    :param _mask_Nx: number of pixels in x-direction [1]
    :param _mask_Ny: number of pixels in y-direction [1]
    :param _grid_nx: number of grids in x-direction
    :param _grid_ny: number of grids in y-direction
    :param _grid_sh: grid shape (0: Circular grids case. 1: Rectangular grids case. 2: 2-D phase grating)
    :param _grid_dx: grid dimension in x-direction, width for rectangular or elliptical grids [m]
    :param _grid_dy: grid dimension in y-direction, height for rectangular or elliptical grids [m]
    :param _grid_angle: tilt angle of the grid [rad]
    :param _mask_x0: horizontal coordinate of the mask [m]
    :param _mask_y0: vertical coordinate of the mask [m]
    :return: transmission (SRWLOptT) type optical element which simulates the PMA
    """

    input_parms = { #MR29092016: Added all input parameters to include in return object:  
        "type": "mask",  
        "refractiveIndex": _delta,  
        "attenuationLength": _atten_len,  
        "maskThickness": _thick,  
        "gridShape": _grid_sh,  
        "horizontalGridDimension": _grid_dx,  
        "verticalGridDimension": _grid_dy,  
        "horizontalGridPitch": _pitch_x,  
         "verticalGridPitch": _pitch_y,  
        "horizontalGridsNumber": _grid_nx,  
        "verticalGridsNumber": _grid_ny,
        "horizontalPixelsNumber": _mask_Nx,  
        "verticalPixelsNumber": _mask_Ny,  
        "gridTiltAngle": _grid_angle,  
        "horizontalSamplingInterval": _hx,  
        "verticalSamplingInterval": _mask_Ny,  
        "horizontalMaskCoordinate": _mask_x0,  
        "verticalMaskCoordinate": _mask_y0,  
    }
    
    # Check if _grid_dy is set by user.
    if _grid_dy == 0:
        _grid_dy = _grid_dx  # An ellipse becomes a circle and a rectangle becomes a square.
    if _grid_sh == 2:
        _grid_dx = _pitch_x  # Change grid_size for 2D grating, grid_size equal to pitch
        _grid_dy = _pitch_y
    # Calculate the range of mask.
    mask_Rx = _hx * _mask_Nx  # mask range in x-direction [m].
    mask_Ry = _hy * _mask_Nx  # mask range in y-direction [m].

    # Calculate the range of grid.
    grid_Rx = _pitch_x * _grid_nx  # grid range in x-direction [m].
    grid_Ry = _pitch_y * _grid_ny  # grid range in y-direction [m].

    # Generate Transmission Optical Element.
    trans_opt = SRWLOptT(_nx=_mask_Nx, _ny=_mask_Ny, _rx=mask_Rx, _ry=mask_Ry, _arTr=None, _extTr=0, _x=0, _y=0)

    # Same data alignment as for wavefront: outer loop vs y, inner loop vs x.
    pointer = 0  # pointer for array trans_opt.arTr
    y = - mask_Ry / 2  # Mask is always centered on the grid, however grid can be shifted.

    #Calculate aux. constants for equations for edges of rectangle.
    xCross1 = None; yCross1 = None #OC27012017
    xCross2 = None; yCross2 = None
    k1 = None; k2 = None; k3 = None; k4 = None
    if(_grid_sh == 1):
        grid_dx_d_sqrt2 = _grid_dx / sqrt(2)
        cosGridAng = cos(_grid_angle)
        sinGridAng = sin(_grid_angle)
        xCross2 = grid_dx_d_sqrt2 * cosGridAng
        #xCross1 = - grid_dx_d_sqrt2 * cosGridAng
        xCross1 = -xCross2
        yCross2 = grid_dx_d_sqrt2 * sinGridAng
        #yCross1 = - grid_dx_d_sqrt2 * sinGridAng
        yCross1 = -yCross2
        k1 = tan(pi / 4 + _grid_angle)
        k2 = -tan(pi / 4 - _grid_angle)
        k4 = tan(pi / 4 + _grid_angle)
        k3 = -tan(pi / 4 - _grid_angle)

        #print(k1, k2, k3, k4)

    for iy in range(_mask_Ny):
        # Calculate the relative position in y.
        # NOTE: Use round to solve the precision issue!
        pitch_num_y = floor(round(y / _pitch_y, 9))
        y_rel = y - (pitch_num_y * _pitch_y) - _mask_y0
        if y_rel >= _pitch_y / 2:
            y_rel -= _pitch_y

        #print(y)
        #print('')
        #print('')

        x = - mask_Rx / 2  # Mask is always centered on the grid, however grid can be shifted.
        for ix in range(_mask_Nx):
            # Calculate the relative position in x.
            # NOTE: Use round to solve the precision issue!
            pitch_num_x = floor(round(x / _pitch_x, 9))
            x_rel = x - (pitch_num_x * _pitch_x) - _mask_x0

            if x_rel >= _pitch_x / 2:
                x_rel -= _pitch_x

            # Initialize the bool parameter.
            inside_hole = False
            phase_shift = False

            # Hartmann hole in an elliptical shape.
            if _grid_sh == 0:
                if (x_rel / _grid_dx) ** 2 + (y_rel / _grid_dy) ** 2 < 1 \
                        and not (round(x_rel - (x - _mask_x0), 9) == 0 and round(y_rel - (y - _mask_y0), 9) == 0) \
                        and abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2:
                    inside_hole = True

            # Hartmann hole in a rectangular shape.
            elif _grid_sh == 1:
                # Calculate the equations for edges of rectangle.
                #OC27012017 (commented-out, calculation of these constants was moved outside of the loop)
                #xCross1 = - _grid_dx / (2 ** 0.5) * math.cos(_grid_angle)
                #yCross1 = - _grid_dx / (2 ** 0.5) * math.sin(_grid_angle)
                #xCross2 = + _grid_dx / (2 ** 0.5) * math.cos(_grid_angle)
                #yCross2 = + _grid_dx / (2 ** 0.5) * math.sin(_grid_angle)
                #k1 = math.tan(math.pi / 4 + _grid_angle)
                #k2 = -math.tan(math.pi / 4 - _grid_angle)
                #k4 = math.tan(math.pi / 4 + _grid_angle)
                #k3 = -math.tan(math.pi / 4 - _grid_angle)

                #print((k2 * x_rel + (yCross2 - k2 * xCross2)), y_rel, (k3 * x_rel + (yCross1 - k3 * xCross1)))
                #print((k1 * x_rel + (yCross1 - k1 * xCross1)), y_rel, (k4 * x_rel + (yCross2 - k4 * xCross2)))
                #print(abs(x - _mask_x0), _pitch_x / 2)
                #print(abs(y - _mask_y0), _pitch_y / 2)
                #print(abs(x), grid_Rx / 2)
                #print(abs(y), grid_Ry / 2)

                if (k2 * x_rel + (yCross2 - k2 * xCross2)) > y_rel > (k3 * x_rel + (yCross1 - k3 * xCross1)) \
                        and (k1 * x_rel + (yCross1 - k1 * xCross1)) > y_rel > (k4 * x_rel + (yCross2 - k4 * xCross2)) \
                        and not (abs(x - _mask_x0) < _pitch_x / 2 and abs(y - _mask_y0) < _pitch_y / 2) \
                        and abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2:
                    inside_hole = True
                    print(y)

            # Grating shearing interferometry in a 2D phase grating.
            elif _grid_sh == 2:
                phase_shift = False
                if (x_rel >= 0 and y_rel < 0) or (x_rel < 0 and y_rel >= 0):
                    phase_shift = True

            else:
                raise ValueError('Unknown shape code.')

            # Give values to trans_opt.arTr.
            if inside_hole and not (_grid_sh == 2):
                trans_opt.arTr[pointer] = 1  # amplitude transmission.  (!) not in physics yet
                trans_opt.arTr[pointer + 1] = 0  # optical path difference. (!) not in physics yet
            else:
                trans_opt.arTr[pointer] = 0  # amplitude transmission.  (!) not in physics yet
                trans_opt.arTr[pointer + 1] = 0  # optical path difference. (!) not in physics yet

            if _grid_sh == 2:
                # Give values to OpT.arTr
                # Make final judgement.
                if phase_shift:
                    trans_opt.arTr[pointer] = exp(-0.5 * _thick / _atten_len)  # amplitude transmission.
                    trans_opt.arTr[pointer + 1] = -_delta * _thick  # optical path difference.
                else:
                    trans_opt.arTr[pointer] = 1  # amplitude transmission.
                    trans_opt.arTr[pointer + 1] = 0  # optical path difference.
                if not (abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2):  # check if it is in grid area
                    trans_opt.arTr[pointer] = 0
                    trans_opt.arTr[pointer + 1] = 0

                    # Shift the pointer by 2.
            pointer += 2

            # Step x by _hx.
            x += _hx

        # Step y by _hy.
        y += _hy

    trans_opt.input_parms = input_parms #MR29092016
    return trans_opt

#****************************************************************************
#Version from L. Huang (moved to srwlib on March 11, 2022)
def srwl_opt_setup_Hartmann_sensor(_delta: float,
#def srwl_opt_setup_mask(_delta: float,
                        _atten_len: float,
                        _thick: float,
                        _hx: float,
                        _hy: float,
                        _pitch_x: float,
                        _pitch_y: float,
                        _mask_nx: int,
                        _mask_ny: int,
                        _grid_nx: int,
                        _grid_ny: int,
                        _grid_sh: int,
                        _grid_dx: float,
                        _grid_dy: float = 0,
                        _grid_angle: float = 0,
                        _mask_x0: float = 0,
                        _mask_y0: float = 0,
                        ):
    """
    Setup Transmission type Optical Element which simulates a mask array for at-wavelength metrology.

    Parameters
    ----------
        _delta: `float`
            The refractive index decrement (can be one number of array vs photon energy)
        _atten_len: `float`
            The attenuation length [m] (can be one number of array vs photon energy)
        _thick: `float`
            The thickness of mask [m]
        _hx: `float`
            The sampling interval in x-direction [m]
        _hy: `float`
            The sampling interval in y-direction [m]
        _pitch_x: `float`
            The grid pitch in x-direction [m]
        _pitch_y: `float`
            The grid pitch in y-direction [m]
        _mask_nx: `int`
            The number of pixels in x-direction [1]
        _mask_ny: `int`
            The number of pixels in y-direction [1]
        _grid_nx: `int`
            The number of grids in x-direction
        _grid_ny: `int`
            The number of grids in y-direction
        _grid_sh: `int`
            The grid shape. 0: Circular grids case. 1: Rectangular grids case. 2: 2-D phase grating
        _grid_dx: `float`
            The grid dimension in x-direction, width for rectangular or elliptical grids [m]
        _grid_dy: `float`
            The grid dimension in y-direction, height for rectangular or elliptical grids [m]
        _grid_angle: `float`
            The tilt angle of the grid [rad]
        _mask_x0: `float`
            The center of mask in x
        _mask_y0: `float`
            The center of mask in y
    Returns
    -------
        `SRWLOptT`
            transmission (SRWLOptT) type optical element which simulates the PMA
    """

    try: #OC13032022
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')

    # Check if _grid_dy is set by user.
    if _grid_dy == 0:
        _grid_dy = _grid_dx  # An ellipse becomes a circle and a rectangle becomes a square.
    if _grid_sh == 2:
        _grid_dx = _pitch_x  #Change grid_size for 2D grating, grid_size eauql to pitch
        _grid_dy = _pitch_y
    # Calculate the range of mask.
    mask_Rx = _hx * (_mask_nx-1)  # mask range in x-direction [m].
    mask_Ry = _hy * (_mask_ny-1)  # mask range in y-direction [m].

    # Calculate the range of grid.
    grid_Rx = _pitch_x * _grid_nx  # grid range in x-direction [m].
    grid_Ry = _pitch_y * _grid_ny  # grid range in y-direction [m].

    if _grid_sh == 1:
        # Calculate the equations for edges of rectangle.
        xTpLt0 = -_grid_dx/2
        yTpLt0 = _grid_dy/2
        xTpLt = xTpLt0 * cos(_grid_angle) - yTpLt0*sin(_grid_angle)
        yTpLt = xTpLt0 * sin(_grid_angle) + yTpLt0*cos(_grid_angle)

        xBmRt0 = _grid_dx/2
        yBmRt0 = -_grid_dy/2
        xBmRt = xBmRt0 * cos(_grid_angle) - yBmRt0*sin(_grid_angle)
        yBmRt = xBmRt0 * sin(_grid_angle) + yBmRt0*cos(_grid_angle)

        kTB = tan(_grid_angle)
        if _grid_angle==0:
            kLR = -1e32
        elif _grid_angle==pi/2:
            kLR = 0.0
        else:
            kLR = -1/kTB

    # Use numpy array
    x, y = np.meshgrid(np.linspace(-mask_Rx/2, mask_Rx/2, num=_mask_nx),
                       np.linspace(-mask_Ry/2, mask_Ry/2, num=_mask_ny))

    # Calculate the relative position in y.
    # NOTE: Use round to solve the precision issue!
    pitch_num_y = np.floor(np.round(y / _pitch_y, 9))
    y_rel = np.round(y - (pitch_num_y * _pitch_y) - _mask_y0, 12)
    y_rel[y_rel>=_pitch_y/2] = np.round(y_rel[y_rel>=_pitch_y/2] - _pitch_y, 12)

    # Calculate the relative position in x.
    # NOTE: Use round to solve the precision issue!
    pitch_num_x = np.floor(np.round(x / _pitch_x, 9))
    x_rel = np.round(x - (pitch_num_x * _pitch_x) - _mask_x0, 12)
    x_rel[x_rel>=_pitch_x/2] = np.round(x_rel[x_rel>=_pitch_x/2] - _pitch_x, 12)

    # Initialize the bool parameter.
    inside_hole = np.zeros(x.shape, dtype=bool)

    # Hartmann subaperture in an elliptical shape.
    if _grid_sh == 0:

        inside_hole[np.logical_and.reduce(
            ((x_rel/(_grid_dx/2))**2 + (y_rel/(_grid_dy/2))**2 <= 1,
             np.logical_not(np.logical_and(
                 np.round(x_rel-(x-_mask_x0), 9) == 0,
                 np.round(y_rel-(y-_mask_y0), 9) == 0)),
                np.abs(x) < grid_Rx / 2,
                np.abs(y) < grid_Ry / 2
             )
        )
        ] = True

    # Hartmann subaperture in a rectangular shape.
    elif _grid_sh == 1:

        inside_hole[np.logical_and.reduce(
            (y_rel < kTB * (x_rel - xTpLt) + yTpLt,
             y_rel > kTB * (x_rel - xBmRt) + yBmRt,
             y_rel < kLR * (x_rel - xBmRt) + yBmRt,
             y_rel > kLR * (x_rel - xTpLt) + yTpLt,
             np.logical_not(np.logical_and(
                 np.abs(x-_mask_x0) < _pitch_x/2,
                 np.abs(y-_mask_y0) < _pitch_y/2)),
             np.abs(x) < grid_Rx / 2,
             np.abs(y) < grid_Ry / 2
             )
        )
        ] = True

    # Grating shearing interferometry in a 2D phase grating.
    elif _grid_sh == 2:

        phase_shift = np.zeros(x.shape, dtype=bool)
        phase_shift[np.logical_or(
            np.logical_and(x_rel >= 0, y_rel < 0),
            np.logical_and(x_rel < 0, y_rel >= 0))] = True

    else:
        print('''Unknown shape code.''')  # (!)

    # Give values to trans_opt.arTr.
    if _grid_sh!=2:
        amp = np.zeros(x.shape) # amplitude transmission.
        amp[inside_hole] = 1 # amplitude transmission.
        opd = np.zeros(x.shape) # optical path difference. (!) not in physics yet

    else:
        amp = np.ones(x.shape) # amplitude transmission.
        amp[phase_shift] = exp(-0.5*_thick/_atten_len) # amplitude transmission.

        opd = np.zeros(x.shape) # optical path difference.
        opd[phase_shift] = -_delta*_thick # optical path difference.

        # If it is outside grid area, ...
        amp[np.logical_not(np.logical_and(np.abs(x) < grid_Rx / 2, np.abs(y) < grid_Ry / 2))] = 0  # optical path difference.
        opd[np.logical_not(np.logical_and(np.abs(x) < grid_Rx / 2, np.abs(y) < grid_Ry / 2))] = 0 # optical path difference.

    arTr_ndarray = np.array([amp.reshape(amp.size), opd.reshape(opd.size)]).transpose()
    #OC13032022
    auxSize = arTr_ndarray.size
    arTr_array = arTr_ndarray.reshape(auxSize)
    #arTr_list = arTr_ndarray.reshape(arTr_ndarray.size).tolist()
    #arTr_array = array('d', arTr_list)

    # Generate Transmission Optical Element.
    trans_opt = SRWLOptT(_nx=_mask_nx, _ny=_mask_ny, _rx=mask_Rx, _ry=mask_Ry, _arTr=arTr_array, _extTr=0, _x=0, _y=0)

    return trans_opt

#****************************************************************************
#The following version is under development 
def srwl_opt_setup_Hartmann_sensor_dev(_per_x, _per_y, _hole_nx, _hole_ny, _hole_sh, _hole_dx, _hole_dy=0, _hole_ang=0, _ang=0, 
#def srwl_opt_setup_Hartmann_sensor(_per_x, _per_y, _hole_nx, _hole_ny, _hole_sh, _hole_dx, _hole_dy=0, _hole_ang=0, _ang=0, 
                                   _delta=0, _atten_len=1e-12, _thick=None, 
                                   _nx=1001, _ny=1001, _rx=0, _ry=0, _e_start=0, _e_fin=0, _xc=0, _yc=0):
    """Setup Transmission type Optical Element which simulates a Harmann sensor, i.e. mask array/matrix for at-wavelength metrology.
    :param _per_x: hole placing period in horiz. direction [m]
    :param _per_y: hole placing period in vert. direction [m]
    :param _hole_nx: number of holes in x-direction
    :param _hole_ny: number of holes in y-direction
    :param _hole_sh: hole shape (0- elliptical, 1- rectangular)
    :param _hole_dx: hole size in horiz. direction (before rotation) [m]
    :param _hole_dy: hole size in vert. direction (before rotation) [m]
    :param _hole_ang: tilt/rotation angle of each hole around its center [rad]
    :param _ang: tilt/rotation angle of entire sensor/mask around (_xc, _yc) [rad]
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _thick: thickness of mask/sensor [m]
    :param _nx: number of points in horiz. direction to represent the transmission object
    :param _ny: number of points in vert. direction to represent the transmission object
    :param _rx: size of the mask/sensor in horiz. direction [m]
    :param _ry: size of the mask/sensor in vert. direction [m]
    :param _e_start: initial photon energy (to be used if _delta and _atten_len are arrays)
    :param _e_fin: final photon energy (to be used if _delta and _atten_len are arrays)
    :param _xc: horizontal coordinate of the mask center [m]
    :param _yc: vertical coordinate of the mask center [m]
    :return: transmission (SRWLOptT) type optical element which simulates the Harmann sensor
    """

    hole_dx = _hole_dx
    hole_dy = _hole_dy
    if(hole_dy == 0): hole_dy = hole_dx

    halfHoleDx = 0.5*hole_dx
    halfHoleDy = 0.5*hole_dy

    rxLoc = _per_x*_hole_nx if(_rx <= 0) else _rx #hor. size of the transmission object
    ryLoc = _per_y*_hole_ny if(_ry <= 0) else _ry #vert. size of the transmission object

    half_rxLoc = 0.5*rxLoc
    half_ryLoc = 0.5*ryLoc

    rx = rxLoc
    ry = ryLoc
    sinAng = 0
    cosAng = 1
    if(_ang != 0):
        sinAng = sin(_ang)
        cosAng = cos(_ang)
        rx = rxLoc*cosAng + ryLoc*sinAng
        ry = rxLoc*sinAng + ryLoc*cosAng

    half_rx = 0.5*rx
    half_ry = 0.5*ry
    xStartGen = _xc - half_rx
    yStartGen = _yc - half_ry

    sinAngHole = 0
    cosAngHole = 1
    if(_hole_ang != 0):
        sinAngHole = sin(_hole_ang)
        cosAngHole = cos(_hole_ang)

    ne = 1
    if((_thick is not None) 
       and (isinstance(_delta, list) or isinstance(_delta, array))
       and (isinstance(_atten_len, list) or isinstance(_atten_len, array))): 
        neDelta = len(_delta)
        neAttLen = len(_atten_len)
        ne = neDelta if(neDelta <= neAttLen) else neAttLen

    delta = None
    attLen = None
    if(ne == 1):
        if(isinstance(_delta, list) or isinstance(_delta, array)): delta = _delta[0]
        if(isinstance(_atten_len, list) or isinstance(_atten_len, array)): attLen = _atten_len[0]

    extTr = 0
    if((_rx > 0) or (_ry > 0)): extTr = 1
    #Transmission Element (initialized according to extTr)
    opTr = SRWLOptT(_nx=_nx, _ny=_ny, _rx=rx, _ry=ry, _extTr=extTr, _x=_xc, _y=_yc, _ne=ne, _eStart=_e_start, _eFin=_e_fin, _alloc_base=[extTr,0])

    xStep = rx/(_nx - 1)
    yStep = ry/(_ny - 1)

    #Determining hole "footprint" dimensions after all rotations
    xSizeHole = hole_dx
    ySizeHole = hole_dy
    angTot = _hole_ang + _ang
    if(abs(angTot) > 1.e-07): #to steer
        sinAngTot = sin(angTot)
        cosAngTot = cos(angTot)
        xSizeHole = hole_dx*cosAngTot + hole_dy*sinAngTot
        ySizeHole = hole_dx*sinAngTot + hole_dy*cosAngTot

    xHalfSizeHole = 0.5*xSizeHole
    yHalfSizeHole = 0.5*ySizeHole

    xHoleCenMinLoc = _xc - 0.5*(_hole_nx - 1)*_per_x
    yHoleCenMinLoc = _yc - 0.5*(_hole_ny - 1)*_per_y

    xPerMesh = 2
    yPerMesh = xPerMesh*_nx

    #Loop over Holes
    yHoleCenLoc = yHoleCenMinLoc
    for ihy in range(_hole_ny):
        xHoleCenLoc = xHoleCenMinLoc
        for ihx in range(_hole_nx):

            #Transform ccordinates of hole center to lab frame, if necessary
            yHoleCen = yHoleCenLoc
            xHoleCen = xHoleCenLoc
            if(_ang != 0):
                xHoleCenLoc_mi_xc = xHoleCenLoc - _xc
                yHoleCenLoc_mi_yc = yHoleCenLoc - _yc
                xHoleCen = _xc + xHoleCenLoc_mi_xc*cosAng - yHoleCenLoc_mi_yc*sinAng #check sign before sinAng and treatment of _xc, _yc
                yHoleCen = _yc + yHoleCenLoc_mi_yc*cosAng + xHoleCenLoc_mi_xc*sinAng

            xStart = xHoleCen - xHalfSizeHole
            dixStart = (xStart - xStartGen)/xStep
            ixStart = int(dixStart)
            if((dixStart - ixStart) > 1e-12): ixStart += 1
            xStartMesh = xStartGen + ixStart*xStep
            xEnd = xHoleCen + xHalfSizeHole
            dixEnd = (xEnd - xStartGen)/xStep
            ixEnd = int(dixEnd)
            if((dixEnd - ixEnd) > 1e-12): ixEnd -= 1
            nxCur = ixEnd - ixStart + 1

            yStart = yHoleCen - yHalfSizeHole
            diyStart = (yStart - yStartGen)/yStep
            iyStart = int(diyStart)
            if((diyStart - iyStart) > 1e-12): iyStart += 1
            yStartMesh = yStartGen + iyStart*yStep
            yEnd = yHoleCen + yHalfSizeHole
            diyEnd = (yEnd - yStartGen)/yStep
            iyEnd = int(diyEnd)
            if((diyEnd - iyEnd) > 1e-12): iyEnd -= 1
            nyCur = iyEnd - iyStart + 1

            #Loop over Transmission points vs x and y around one hole
            y = yStartMesh
            for iy in range(nyCur):
                x = xStartMesh
                for ix in range(nxCur):

                    #Transform (x, y) to the Mask frame
                    xLocMask = x
                    yLocMask = y
                    if(_ang != 0):
                        x_mi_xc = x - _xc
                        y_mi_yc = y - _yc
                        xLocMask = _xc + x_mi_xc*cosAng + y_mi_yc*sinAng #check sign before sinAng and treatment of _xc, _yc
                        yLocMask = _yc + y_mi_yc*cosAng - x_mi_xc*sinAng

                    #Transform (xLocMask, yLocMask) to the Hole frame
                    xLocHole = xLocMask
                    yLocHole = yLocMask
                    if(_hole_ang != 0):
                        xLocMask_mi_xHoleCen = xLocMask - xHoleCen
                        yLocMask_mi_yHoleCen = yLocMask - yHoleCen
                        xLocHole = xHoleCen + xLocMask_mi_xHoleCen*cosAngHole + yLocMask_mi_yHoleCen*sinAngHole #check sign before sin and treatment of center point coords
                        yLocHole = yHoleCen + yLocMask_mi_yHoleCen*cosAngHole - xLocMask_mi_xHoleCen*sinAngHole
                    
                    xLocHole_mi_xHoleCen = xLocHole - xHoleCenLoc #check if this is necessary
                    yLocHole_mi_yHoleCen = yLocHole - yHoleCenLoc #check if this is necessary

                    otst = -1
                    if(_hole_sh == 0): #Elliptical hole shape
                        relX = xLocHole_mi_xHoleCen/halfHoleDx
                        relY = yLocHole_mi_yHoleCen/halfHoleDy
                        if((relX*relX + relY*relY) < 1): ofst = iy*yPerMesh + ix*xPerMesh
                    elif(_hole_sh == 1): #Rectangular hole shape
                        if((-halfHoleDx < xLocHole_mi_xHoleCen) and (xLocHole_mi_xHoleCen < halfHoleDx) and
                           (-halfHoleDy < yLocHole_mi_yHoleCen) and (yLocHole_mi_yHoleCen < halfHoleDy)): 
                            ofst = iy*yPerMesh + ix*xPerMesh

                    if(otst >= 0):
                        if(ne == 1):
                            opTr.arTr[ofst] = 1.
                            opTr.arTr[ofst + 1] = 0.
                        #elif(ne > 1):
                        #ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd


                    x += xStep
                y += yStep

            xHoleCenLoc += _per_x
        yHoleCenLoc += _per_y

    return opTr

#****************************************************************************
def srwl_opt_setup_surf_height_1d(_height_prof_data, _dim, _ang, _ang_r=0, _amp_coef=1, _ar_arg_long=None, _nx=0, _ny=0, _size_x=0, _size_y=0, _xc=0, _yc=0): #OC19072018
    """
    Setup Transmission type optical element with 1D (mirror or grating) surface Heght Profile data
    :param _height_prof_data: two- or one-column table containing, in case of two columns: longitudinal position in [m] (1st column) and the Height Profile in [m] (2nd column) data; in case of one column, it contains the Height Profile data
    :param _dim: orientation of the reflection (deflection) plane; can be 'x' or 'y'
    :param _ang: grazing angle (between input optical axis and mirror/grating plane)
    :param _ang_r: reflection angle (between output optical axis and mirror/grating plane)
    :param _amp_coef: height profile "amplification coefficient"
    :param _ar_arg_long: optional array of longitudinal position (along mirror/grating) in [m]; if _ar_arg_long is not None, any longitudinal position contained in _height_prof_data is ignored
    :param _nx: optional number of points in horizontal dimension of the output transmission optical element
    :param _ny: optional number of points in vertical dimension of the output transmission optical element
    :param _size_x: optional horizontal transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :param _size_y: optional vertical transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :param _xc: optional horizontal center position of the output transmission optical element
    :param _yc: optional vertical center position of the output transmission optical element
    :return: transmission (SRWLOptT) type optical element which simulates the effect of surface height error
    """
    #To test all options!

    input_parms = { #MR26022016: Added all input parameters to include in return object:
        "type": "mirror",
        "heightProfileFile": "",
        "orientation": _dim,
        "grazingAngle": _ang,
        "reflectionAngle": _ang_r,
        "heightAmplification": _amp_coef,
        "longitudinalPosition": _ar_arg_long,
        "horizontalPoints": _nx,
        "verticalPoints": _ny,
        "horizontalTransverseSize": _size_x,
        "verticalTransverseSize": _size_y,
        "horizontalCenterPosition": _xc,
        "verticalCenterPosition": _yc,
    }

    if((_dim != 'x') and (_dim != 'y')): raise Exception('The orientation of the deflection plane can be either horizontal (\'x\') or vertical (\'y\')')

    if(_ang_r == 0): _ang_r = _ang
    sinAng = sin(_ang)
    sinAngR = sin(_ang_r)

    if _ar_arg_long is None:
        argHeightProfData = _height_prof_data[0]
        valHeightProfData = _height_prof_data[1]
    else:
        argHeightProfData = _ar_arg_long
        if len(_height_prof_data) >= 2:
            valHeightProfData = _height_prof_data[1]
        else: valHeightProfData = _height_prof_data[0]

    npData = len(valHeightProfData)
    npDataTr = 100 #default value
    sizeLongProj = (argHeightProfData[npData - 1] - argHeightProfData[0])*sinAngR
    sizeTr = 1 #default value

    nx = _nx
    if nx <= 0:
        if('x' in _dim): nx = npData
        else: nx = npDataTr
    ny = _ny
    if ny <= 0:
        if('y' in _dim): ny = npData
        else: ny = npDataTr

    if _size_x > 0: sizeX = _size_x
    else:
        sizeX = sizeLongProj
        if('y' in _dim): sizeX = sizeTr

    if _size_y > 0: sizeY = _size_y
    else:
        sizeY = sizeTr
        if('y' in _dim): sizeY = sizeLongProj

    optSlopeErr = SRWLOptT(nx, ny, sizeX, sizeY, _x=_xc, _y=_yc)

    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)

    perX = 2
    perY = 2*auxMesh.nx

    hApprox = 0
    ipStart = 0

    #DEBUG
    #print('Heigh profile: nx=', auxMesh.nx, ' ny=', auxMesh.ny)

    #if('y' in _dim):
    arg1Start = auxMesh.yStart - _yc
    arg2Start = auxMesh.xStart - _xc
    n1 = auxMesh.ny; n2 = auxMesh.nx
    per1 = perY; per2 = perX
    arg1Step = yStep; arg2Step = xStep
    
    if('x' in _dim):
        arg1Start = auxMesh.xStart - _xc
        arg2Start = auxMesh.yStart - _yc
        n1 = auxMesh.nx; n2 = auxMesh.ny
        per1 = perX; per2 = perY
        arg1Step = xStep; arg2Step = yStep

    arg1 = arg1Start #to make sure that only the mesh moves
    for i1 in range(n1):
        hApprox = 0
        arg1_1 = argHeightProfData[ipStart]*sinAngR
        for i in range(ipStart + 1, npData):
            arg1_2 = argHeightProfData[i]*sinAngR
            if((arg1_1 <= arg1) and (arg1 < arg1_2)):
                hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAngR))*(arg1 - arg1_1) + valHeightProfData[i-1]
                ipStart = i - 1
                break
            arg1_1 = arg1_2

        per1_i1 = per1*i1
        arg2 = arg2Start #to make sure that only the mesh moves
        for i2 in range(n2):
            ofst = per2*i2 + per1_i1
            optSlopeErr.arTr[ofst] = 1. #Amplitude Transmission #consider taking into account reflectivity
            optSlopeErr.arTr[ofst + 1] = 0. #Optical Path Difference
            if(hApprox != 0):
                optSlopeErr.arTr[ofst + 1] = -(sinAng + sinAngR)*hApprox*_amp_coef #Optical Path Difference (to check sign!)

            arg2 += arg2Step
        arg1 += arg1Step

    optSlopeErr.input_parms = input_parms  #MR16022016
    return optSlopeErr

#****************************************************************************
def srwl_opt_setup_surf_height_1d_old(_height_prof_data, _dim, _ang, _ang_r=0, _amp_coef=1, _ar_arg_long=None, _nx=0, _ny=0, _size_x=0, _size_y=0, _xc=0, _yc=0):
    #OC19072018: This version is very slow for _dim = 'x' with Py 2.7
    """
    Setup Transmission type optical element with 1D (mirror or grating) surface Heght Profile data
    :param _height_prof_data: two- or one-column table containing, in case of two columns: longitudinal position in [m] (1st column) and the Height Profile in [m] (2nd column) data; in case of one column, it contains the Height Profile data
    :param _dim: orientation of the reflection (deflection) plane; can be 'x' or 'y'
    :param _ang: grazing angle (between input optical axis and mirror/grating plane)
    :param _ang_r: reflection angle (between output optical axis and mirror/grating plane)
    :param _amp_coef: height profile "amplification coefficient"
    :param _ar_arg_long: optional array of longitudinal position (along mirror/grating) in [m]; if _ar_arg_long is not None, any longitudinal position contained in _height_prof_data is ignored
    :param _nx: optional number of points in horizontal dimension of the output transmission optical element
    :param _ny: optional number of points in vertical dimension of the output transmission optical element
    :param _size_x: optional horizontal transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :param _size_y: optional vertical transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :param _xc: optional horizontal center position of the output transmission optical element
    :param _yc: optional vertical center position of the output transmission optical element
    :return: transmission (SRWLOptT) type optical element which simulates the effect of surface height error
    """
    #To test all options!

    input_parms = { #MR26022016: Added all input parameters to include in return object:
        "type": "mirror",
        "heightProfileFile": "",
        "orientation": _dim,
        "grazingAngle": _ang,
        "reflectionAngle": _ang_r,
        "heightAmplification": _amp_coef,
        "longitudinalPosition": _ar_arg_long,
        "horizontalPoints": _nx,
        "verticalPoints": _ny,
        "horizontalTransverseSize": _size_x,
        "verticalTransverseSize": _size_y,
        "horizontalCenterPosition": _xc,
        "verticalCenterPosition": _yc,
    }

    if(_ang_r == 0): _ang_r = _ang
    sinAng = sin(_ang)
    sinAngR = sin(_ang_r)

    if _ar_arg_long is None:
        argHeightProfData = _height_prof_data[0]
        valHeightProfData = _height_prof_data[1]
    else:
        argHeightProfData = _ar_arg_long
        if len(_height_prof_data) >= 2:
            valHeightProfData = _height_prof_data[1]
        else: valHeightProfData = _height_prof_data[0]

    npData = len(valHeightProfData)
    npDataTr = 100 #default value
    sizeLongProj = (argHeightProfData[npData - 1] - argHeightProfData[0])*sinAngR
    sizeTr = 1 #default value

    nx = _nx
    if nx <= 0:
        if('x' in _dim): nx = npData
        else: nx = npDataTr
    ny = _ny
    if ny <= 0:
        if('y' in _dim): ny = npData
        else: ny = npDataTr

    if _size_x > 0: sizeX = _size_x
    else:
        sizeX = sizeLongProj
        if('y' in _dim): sizeX = sizeTr

    if _size_y > 0: sizeY = _size_y
    else:
        sizeY = sizeTr
        if('y' in _dim): sizeY = sizeLongProj

    #optSlopeErr = SRWLOptT(nx, ny, sizeX, sizeY)
    optSlopeErr = SRWLOptT(nx, ny, sizeX, sizeY, _x=_xc, _y=_yc)

    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)

    #y = auxMesh.yStart
    y = auxMesh.yStart - _yc #to make sure that only the mesh moves
    hApprox = 0
    ipStart = 0

    #DEBUG
    #print('Heigh profile: nx=', auxMesh.nx, ' ny=', auxMesh.ny)
    
    #for iy in range(optSlopeErr.ny):
    for iy in range(auxMesh.ny):
        if('y' in _dim):
            hApprox = 0
            y1 = argHeightProfData[ipStart]*sinAngR
            for i in range(ipStart + 1, npData):
                y2 = argHeightProfData[i]*sinAngR
                if((y1 <= y) and (y < y2)):
                    hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAngR))*(y - y1) + valHeightProfData[i-1]
                    #hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAng))*(y - y1) + valHeightProfData[i-1]
                    #print(ipStart, i, iy, y1, y, y2, argHeightProfData[i-1], argHeightProfData[i], valHeightProfData[i-1], valHeightProfData[i], hApprox)
                    ipStart = i - 1
                    break
                y1 = y2

        #x = auxMesh.xStart
        x = auxMesh.xStart - _xc #to make sure that only the mesh moves
        
        #for ix in range(optSlopeErr.nx):
        for ix in range(auxMesh.nx):
            if('x' in _dim):
                if(ix == 0): ipStart = 0
                hApprox = 0
                x1 = argHeightProfData[ipStart]*sinAngR
                for i in range(ipStart + 1, npData):
                    x2 = argHeightProfData[i]*sinAngR
                    if((x1 <= x) and (x < x2)):
                        hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAngR))*(x - x1) + valHeightProfData[i-1]
                        #hApprox = ((valHeightProfData[i] - valHeightProfData[i-1])/((argHeightProfData[i] - argHeightProfData[i-1])*sinAng))*(x - x1) + valHeightProfData[i-1]
                        ipStart = i - 1
                        break
                    x1 = x2
            #ofst = 2*ix + (2*optSlopeErr.nx)*iy
            ofst = 2*ix + (2*auxMesh.nx)*iy

            optSlopeErr.arTr[ofst] = 1. #Amplitude Transmission #consider taking into account reflectivity
            optSlopeErr.arTr[ofst + 1] = 0. #Optical Path Difference
            if(hApprox != 0):
                #optSlopeErr.arTr[ofst + 1] = -2*sinAng*hApprox #Optical Path Difference (to check sign!)
                optSlopeErr.arTr[ofst + 1] = -(sinAng + sinAngR)*hApprox*_amp_coef #Optical Path Difference (to check sign!)
                #print(ix, iy, optSlopeErr.arTr[ofst + 1])
            x += xStep
        y += yStep

    optSlopeErr.input_parms = input_parms  #MR16022016
    return optSlopeErr

#****************************************************************************
def srwl_opt_setup_surf_height_2d(_height_prof_data, _dim, _ang, _ang_r=0, _amp_coef=1, _ar_arg_long=None, _ar_arg_tr=None, _nx=0, _ny=0, _size_x=0, _size_y=0):
    """
    Setup Transmission type optical element with 2D (mirror or grating) surface Heght Profile data
    :param _height_prof_data: a matrix (2D array) containing the Height Profile data in [m]; if _ar_height_prof_x is None and _ar_height_prof_y is None: the first column in _height_prof_data is assumed to be the "longitudinal" position [m] and first row the "transverse" position [m], and _height_prof_data[0][0] is not used; otherwise the "longitudinal" and "transverse" positions on the surface are assumed to be given by _ar_height_prof_x, _ar_height_prof_y 
    :param _dim: orientation of the reflection (deflection) plane; can be 'x' or 'y'
    :param _ang: grazing angle (between input optical axis and mirror/grating plane)
    :param _ang_r: reflection angle (between output optical axis and mirror/grating plane)
    :param _amp_coef: height profile "amplification coefficient"
    :param _ar_arg_long: optional array of longitudinal position (along mirror/grating) in [m] 
    :param _ar_arg_tr: optional array of transverse position on mirror/grating surface in [m] 
    :param _nx: optional number of points in horizontal dimension of the output transmission optical element
    :param _ny: optional number of points in vertical dimension of the output transmission optical element
    :param _size_x: optional horizontal transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :param _size_y: optional vertical transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :return: transmission (SRWLOptT) type optical element which simulates the effect of surface height error
    """
    #To test all options!

    input_parms = { #MR26022016: Options will be used for 2D mirror profiles in Sirepo in the future:
        #"type": "mirror",
        "type": "mirror2d",
        "heightProfileFile": "",
        "orientation": _dim,
        "grazingAngle": _ang,
        "reflectionAngle": _ang_r,
        "heightAmplification": _amp_coef,
        "longitudinalPosition": _ar_arg_long,
        "transversePosition": _ar_arg_tr,
        "horizontalPoints": _nx,
        "verticalPoints": _ny,
        "horizontalTransverseSize": _size_x,
        "verticalTransverseSize": _size_y,
    }

    if(_ang_r == 0): _ang_r = _ang
    sinAng = sin(_ang)
    sinAngR = sin(_ang_r)
    
    #argHeightProfData = _ar_arg_long
    if _ar_arg_long is None:
        npData = len(_height_prof_data[0]) - 1
        sizeLong = _height_prof_data[0][npData] - _height_prof_data[0][1] #OC17122018
        #sizeLong = _height_prof_data[0][npData - 1] - _height_prof_data[0][1]
    else:
        npData = len(_ar_arg_long)
        sizeLong = _ar_arg_long[npData - 1] - _ar_arg_long[0]
        
    sizeLongProj = sizeLong*sinAngR

    if _ar_arg_tr is None:
        npDataTr = len(_height_prof_data) - 1
        sizeTr = _height_prof_data[npDataTr][0] - _height_prof_data[1][0] #OC17122018
        #sizeTr = _height_prof_data[npDataTr - 1][0] - _height_prof_data[1][0]
    else:
        npDataTr = len(_ar_arg_tr)
        sizeTr = _ar_arg_tr[npDataTr - 1] - _ar_arg_tr[0]

    #npData = len(_height_prof_data[0])
    #npDataTr = len(_height_prof_data)

    nx = _nx
    if nx <= 0:
        if('x' in _dim): nx = npData
        else: nx = npDataTr
    ny = _ny
    if ny <= 0:
        if('y' in _dim): ny = npData
        else: ny = npDataTr

    if _size_x > 0: sizeX = _size_x
    else:
        sizeX = sizeLongProj
        if('y' in _dim): sizeX = sizeTr

    if _size_y > 0: sizeY = _size_y
    else:
        sizeY = sizeTr
        if('y' in _dim): sizeY = sizeLongProj

    #sizeX = sizeLongProj; sizeY = sizeTr
    #if('y' in _dim):
    #    sizeX = sizeTr; sizeY = sizeLongProj

    #print(npData, npDataTr)
    #print(sizeLong, sizeTr)
    #print(sizeLongProj, sizeTr)
    #print(sizeX, sizeY)

    optSlopeErr = SRWLOptT(nx, ny, sizeX, sizeY)
    
    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)

    #print(auxMesh.xStart, auxMesh.xFin, auxMesh.nx, xStep)
    #print(auxMesh.yStart, auxMesh.yFin, auxMesh.ny, yStep)

    xTolEdge = 0.001*abs(xStep) #OC18122018
    yTolEdge = 0.001*abs(yStep)

    y = auxMesh.yStart
    hApprox = 0
    
    ipStart = 1
    ipStartTr = 1

    for iy in range(auxMesh.ny):
        y1 = 0; y2 = 0

        if('y' in _dim):
            ipStartTr = 1
            #y1 = argHeightProfData[ipStart]*sinAngR
            if _ar_arg_long is None: y1 = _height_prof_data[0][ipStart]*sinAngR
            else: y1 = _ar_arg_long[ipStart - 1]*sinAngR
            
            for i in range(ipStart + 1, npData + 1):
            #for i in range(ipStart + 1, npData):
                #y2 = argHeightProfData[i]*sinAngR
                if _ar_arg_long is None: y2 = _height_prof_data[0][i]*sinAngR
                else: y2 = _ar_arg_long[i - 1]*sinAngR
                
                if((y1 <= y) and (y < y2)):
                    ipStart = i - 1
                    break
                y1 = y2

            if(y1 == y2): #OC18122018
                if((abs(y - y2) < yTolEdge) and (abs(y2 - auxMesh.yFin) < yTolEdge)):
                    y1 = auxMesh.yFin - yStep

        elif('x' in _dim):
            ipStart = 1
            if _ar_arg_tr is None: y1 = _height_prof_data[ipStartTr][0]
            else: y1 = _ar_arg_tr[ipStartTr - 1]
            
            for i in range(ipStartTr + 1, npDataTr + 1):
            #for i in range(ipStartTr + 1, npDataTr):
                if _ar_arg_tr is None: y2 = _height_prof_data[i][0]
                else: y2 = _ar_arg_tr[i - 1]
                
                if((y1 <= y) and (y < y2)):
                    ipStartTr = i - 1
                    break
                y1 = y2

            if(y1 == y2): #OC18122018
                if((abs(y - y2) < yTolEdge) and (abs(y2 - auxMesh.yFin) < yTolEdge)):
                    y1 = auxMesh.yFin - yStep

        x = auxMesh.xStart
        for ix in range(auxMesh.nx):
            x1 = 0; x2 = 0

            if('y' in _dim):
                if(ix == 0): ipStartTr = 1
              
                if _ar_arg_tr is None: x1 = _height_prof_data[ipStartTr][0]
                else: x1 = _ar_arg_tr[ipStartTr - 1]

                #print(ipStartTr + 1, npDataTr + 1)
                
                for i in range(ipStartTr + 1, npDataTr + 1):
                #for i in range(ipStartTr + 1, npDataTr):
                    if _ar_arg_tr is None: x2 = _height_prof_data[i][0]
                    else: x2 = _ar_arg_tr[i - 1]
                    
                    if((x1 <= x) and (x < x2)):
                        ipStartTr = i - 1
                        #print(ix, iy, x1, x, x2)
                        break
                    #print(ix, i, x1, x2, x)
                    x1 = x2

                if(x1 == x2): #OC18122018
                    if((abs(x - x2) < xTolEdge) and (abs(x2 - auxMesh.xFin) < xTolEdge)):
                        x1 = auxMesh.xFin - xStep
                    
            elif('x' in _dim):
                if(ix == 0): ipStart = 1
                
                #x1 = argHeightProfData[ipStart]*sinAngR
                if _ar_arg_long is None: x1 = _height_prof_data[0][ipStart]*sinAngR
                else: x1 = _ar_arg_long[ipStart - 1]*sinAngR
                
                for i in range(ipStart + 1, npData + 1):
                #for i in range(ipStart + 1, npData):
                    #x2 = argHeightProfData[i]*sinAngR
                    if _ar_arg_long is None: x2 = _height_prof_data[0][i]*sinAngR
                    else: x2 = _ar_arg_long[i - 1]*sinAngR
                    
                    if((x1 <= x) and (x < x2)):
                        ipStart = i - 1
                        break
                    x1 = x2

                if(x1 == x2): #OC18122018
                    if((abs(x - x2) < xTolEdge) and (abs(x2 - auxMesh.xFin) < xTolEdge)):
                        x1 = auxMesh.xFin - xStep

            if _ar_arg_long is not None: ipStart -= 1
            if _ar_arg_tr is not None: ipStartTr -= 1

            #Bi-Linear Interpolation
            xt = 0; yt = 0
            f10 = 0; f01 = 0; f11 = 0
            if(x2 != x1):
                xt = (x - x1)/(x2 - x1)
                if('x' in _dim):
                    f10 = _height_prof_data[ipStartTr][ipStart+1]
                else:
                    f10 = _height_prof_data[ipStartTr+1][ipStart]
            if(y2 != y1):
                yt = (y - y1)/(y2 - y1)
                if('y' in _dim):
                    f01 = _height_prof_data[ipStartTr][ipStart+1]
                else:
                    #OC_TEST
                    #print('ipStartTr=', ipStartTr, 'ipStart=', ipStart)
                    f01 = _height_prof_data[ipStartTr+1][ipStart]
            if((x2 != x1) and (y2 != y1)): f11 = _height_prof_data[ipStartTr+1][ipStart+1]
            
            f00 = _height_prof_data[ipStartTr][ipStart]
            #f10 = heightProfData[ipStartTr+1][ipStart]
            #f01 = heightProfData[ipStartTr][ipStart+1]
            #f11 = heightProfData[ipStartTr+1][ipStart+1]
            a01 = f01 - f00
            a10 = f10 - f00
            a11 = f00 - f01 - f10 + f11
            hApprox = xt*(a10 + a11*yt) + a01*yt + f00

            #print('     x:', x1, x, x2, 'y:', y1, y, y2)
            #print('h:', hApprox, f00, f10, f01, f11, 'ii:', ipStartTr, ipStart)
            #print(' ')

            ofst = 2*ix + (2*auxMesh.nx)*iy
            optSlopeErr.arTr[ofst] = 1. #Amplitude Transmission
            optSlopeErr.arTr[ofst + 1] = 0. #Optical Path Difference
            if(hApprox != 0):
                #optSlopeErr.arTr[ofst + 1] = -2*sinAng*hApprox #Optical Path Difference (to check sign!)
                optSlopeErr.arTr[ofst + 1] = -(sinAng + sinAngR)*hApprox*_amp_coef #Optical Path Difference (to check sign!)
                #print(ix, iy, optSlopeErr.arTr[ofst + 1])
            x += xStep
        y += yStep

    optSlopeErr.input_parms = input_parms #MR26022016
    return optSlopeErr

#****************************************************************************
def srwl_opt_setup_bumps(_ampl, _sx, _sy, _n, _delta, _atten_len, _rx, _ry, _xc=0, _yc=0, _nx=1001, _ny=1001, _n_sig=4, _ampl_min=None, _sx_min=None, _sy_min=None, _seed=None):
    """
    Setup Transmission type Optical Element which simulates a set of Gaussian-shape "Bumps" randomly placed in transverse plane
    :param _ampl: amplitude of bumps [m] or list of min. and max. amplitudes
    :param _sx: horizontal FWHM size of bumps [m] or list of min. and max. horizontal sizes
    :param _sy: vertical FWHM size of bumps [m] or list of min. and max. vertical sizes
    :param _n: number of bumps or list of numbers of bumps of different types
    :param _delta: refractive index decrement of bump material or list of refractive index decrements of bump materials
    :param _atten_len: attenuation length of bump material [m] or list of attenuation lengths of bump materials
    :param _rx: horizontal coordiate range over which bumps are distributed [m]
    :param _ry: vertical coordiate range over which bumps are distributed [m]
    :param _xc: horizontal coordinate of center [m] of the interval / range over bumps are applied [m]
    :param _yc: vertical coordinate of center [m] of the interval / range over bumps are applied [m]
    :param _nx: number of points vs horizontal position to represent the transmission element
    :param _ny: number of points vs vertical position to represent the transmission element
    :param _n_sig: number of sigmas of each Gaussian to take into acount on each side of center position (i.e. only these x, y values will be used:  -_n_sig*sigX < x < _n_sig*sigX, -_n_sig*sigY < y < _n_sig*sigY)
    :param _ampl_min: minimal amplitude of bumps [m]; if it defined, _ampl is treated as maximal amplitude, and the bump amplitude is evenly and randomly distributed between the two values
    :param _sx_min: minimal horizontal FWHM size of bumps [m]; if it defined, _sx is treated as maximal horizontal size, and the bump horizontal size is evenly and randomly distributed between the two values
    :param _sy_min: minimal vertical FWHM size of bumps [m]; if it defined, _sy is treated as maximal vertical size, and the bump vertical size is evenly and randomly distributed between the two values
    :param _seed: integer number to be used to seed random placing of bumps (in None, the seeding will be made from current time)
    :return: transmission (SRWLOptT) type optical element which simulates a set of Gaussian-shape "Bumps" randomly placed in transverse plane
    """

    def SortPair(_pair, _mult=1):
        x1 = _pair[0]*_mult
        x2 = _pair[1]*_mult
        if(x1 > x2):
            aux = x1
            x1 = x2
            x2 = aux
        return [x1, x2]

    multRMS = 1./(2.*sqrt(2.*log(2.)))

    amplIsList = (isinstance(_ampl, list) or isinstance(_ampl, array))
    sizeIsList = ((isinstance(_sx, list) or isinstance(_sx, array)) and (isinstance(_sy, list) or isinstance(_sy, array)))
    nIsList = (isinstance(_n, list) or isinstance(_n, array))
    deltaIsList = (isinstance(_delta, list) or isinstance(_delta, array))
    attLenIsList = (isinstance(_atten_len, list) or isinstance(_atten_len, array))

    if((amplIsList or sizeIsList or deltaIsList or attLenIsList) and (not nIsList)):
        raise Exception("Inconsistent definition of numbers of bumps of different type (_n may need to be a list / array)") 

    arN = _n if(nIsList) else [_n]
    lenArN = len(arN)

    arAmpl = _ampl 
    if(not amplIsList):
        arAmpl = [[_ampl,_ampl]]*lenArN if(_ampl_min is None) else [[_ampl_min,_ampl]]*lenArN
        #TEST
        #print(arAmpl)

    arSx = _sx 
    arSy = _sy 
    if(not sizeIsList): 
        arSx = [[_sx,_sx]]*lenArN if(_sx_min is None) else [[_sx_min,_sx]]*lenArN
        arSy = [[_sy,_sy]]*lenArN if(_sy_min is None) else [[_sy_min,_sy]]*lenArN

    #TEST
    #print('arSx=', arSx)

    arAmplAux = []
    arSxAux = []
    arSyAux = []
    for i in range(lenArN):
        arAmplAux.append(SortPair(arAmpl[i]))
        arSxAux.append(SortPair(arSx[i], multRMS))
        arSyAux.append(SortPair(arSy[i], multRMS))

    arAmpl = arAmplAux
    arSx = arSxAux
    arSy = arSyAux

    #TEST
    #print('arSx', arSx)

    arDelta = _delta if(deltaIsList) else [_delta]*lenArN
    arAttLen = _atten_len if(attLenIsList) else [_atten_len]*lenArN

    if(len(arAmpl) != lenArN): raise Exception("Inconsistent definition of bump amplitudes") 
    if((len(arSx) != lenArN) or (len(arSy) != lenArN)): raise Exception("Inconsistent definition of bump sizes") 
    if(len(arDelta) != lenArN): raise Exception("Inconsistent definition of bump material refractive index decrement(s)") 
    if(len(arAttLen) != lenArN): raise Exception("Inconsistent definition of bump material attenuation length(s)") 

    #sigXmax = abs(_sx)*multRMS
    #sigXmin = sigXmax
    #if(_sx_min is not None): sigXmin = abs(_sx_min)*multRMS
    #if(sigXmin > sigXmax):
    #    aux = sigXmin
    #    sigXmin = sigXmax
    #    sigXmax = aux

    #sigYmax = abs(_sy)*multRMS
    #sigYmin = sigYmax
    #if(_sy_min is not None): sigYmin = abs(_sy_min)*multRMS
    #if(sigYmin > sigYmax):
    #    aux = sigYmin
    #    sigYmin = sigYmax
    #    sigYmax = aux

    #aMax = abs(_ampl)
    #aMin = aMax
    #if(_ampl_min is not None): aMin = abs(_ampl_min)
    #if(aMin > aMax):
    #    aux = aMin
    #    aMin = aMax
    #    aMax = aux

    halfRx = 0.5*_rx
    xStart = _xc - halfRx
    xEnd = _xc + halfRx
    xStep = _rx/(_nx - 1)
    halfRy = 0.5*_ry
    yStart = _yc - halfRy
    yEnd = _yc + halfRy
    yStep = _ry/(_ny - 1)

    perX = 2
    perY = perX*_nx
    op = SRWLOptT(_nx, _ny, _rx, _ry, _alloc_base=[1,0])
    #op = SRWLOptT(_nx, _ny, _rx, _ry)
    #DEBUG
    #print('srwl_opt_setup_bumps: op.arTr memory allocation finished, len(op.arTr)=', len(op.arTr))
    #for i in range(20): print(op.arTr[i])

    if(_seed is None): random.seed(datetime.datetime.now())
    else: random.seed(_seed)

    #ofst = 0
    #for iy in range(_ny):
    #    for ix in range(_nx):
    #        op.arTr[ofst] = 1 #amplitude transmission
    #        op.arTr[ofst + 1] = 0 #optical path difference
    #        ofst += 2
    #DEBUG
    #print('srwl_opt_setup_bumps: initialization of op.arTr finished')

    for ii in range(lenArN):
        #invAttenLen = 1./_atten_len
        invAttenLen = 1./arAttLen[ii]
        delta = arDelta[ii]

        amplMinMax = arAmpl[ii]
        aMin = amplMinMax[0]
        aMax = amplMinMax[1]

        sigXMinMax = arSx[ii]
        sigXmin = sigXMinMax[0]
        sigXmax = sigXMinMax[1]

        sigYMinMax = arSy[ii]
        sigYmin = sigYMinMax[0]
        sigYmax = sigYMinMax[1]

        nBumps = arN[ii]

        #DEBUG
        #print('aMin=', aMin, ' aMax=', aMax)
        #print('sigXmin=', sigXmin, ' sigXmax=', sigXmax)
        #print('sigYmin=', sigYmin, ' sigYmax=', sigYmax)
        #print(nBumps)
        #print('delta=', delta, ' AttLen=', arAttLen[ii])

        for i in range(nBumps):
        #for i in range(_n):
            a = aMin if(aMin == aMax) else random.uniform(aMin, aMax)
            xb = random.uniform(xStart, xEnd)
            yb = random.uniform(yStart, yEnd)
            xbSig = sigXmin if(sigXmin == sigXmax) else random.uniform(sigXmin, sigXmax)
            ybSig = sigYmin if(sigYmin == sigYmax) else random.uniform(sigYmin, sigYmax)
            xbHalfR = _n_sig*xbSig
            ybHalfR = _n_sig*ybSig

            inv2XbSigE2 = 0.5/(xbSig*xbSig)
            inv2YbSigE2 = 0.5/(ybSig*ybSig)

            ixStart = int((xb - xbHalfR - xStart)/xStep + 1.e-09)
            if(ixStart < 0): ixStart = 0
            ixEnd = int((xb + xbHalfR - xStart)/xStep + 1.e-09) + 1
            if(ixEnd >= _nx): ixEnd = _nx

            iyStart = int((yb - ybHalfR - yStart)/yStep + 1.e-09)
            if(iyStart < 0): iyStart = 0
            iyEnd = int((yb + ybHalfR - yStart)/yStep + 1.e-09) + 1
            if(iyEnd >= _ny): iyEnd = _ny

            dy = yStart + iyStart*yStep - yb
            for iy in range(iyStart, iyEnd):

                argY = -dy*dy*inv2YbSigE2
                iy_perY = iy*perY

                dx = xStart + ixStart*xStep - xb
                for ix in range(ixStart, ixEnd):

                    dPath = a*exp(argY - dx*dx*inv2XbSigE2)
                    ofst = iy_perY + ix*perX
                    op.arTr[ofst] *= exp(-0.5*dPath*invAttenLen) #amplitude transmission
                    op.arTr[ofst + 1] += -delta*dPath #optical path difference

                    dx += xStep
                dy += yStep
    return op

#****************************************************************************
def srwl_opt_setup_transit_reg(_delta1, _atten_len1, _delta2, _atten_len2, _thick, _rx, _ry, _w=0, _tr_typ=1, _x0=0, _y0=0, _ang=0, _xc=0, _yc=0, _nx=1001, _ny=1001): #, _e_start=0, _e_fin=0):
    """
    Setup Transmission type Optical Element simulating transition region cutting transverse plane into two half-planes with constant refraction / absorption properties
    :param _delta1: refractive index decrement of one part of the area
    :param _atten_len1: attenuation length of one part of the area
    :param _delta2: refractive index decrement of the other part of the area
    :param _atten_len2: attenuation length of the other part of the area
    :param _thick: thickness of the transmission element [m]
    :param _rx: horizontal coordiate range of the transmission element [m]
    :param _ry: vertical coordiate range of the transmission element [m]
    :param _w: width of transition region between two areas [m]
    :param _tr_typ: type of the transition: 1- linear, 2- ...
    :param _x0: horizontal coordinate of a point on the splitting / transition line [m]
    :param _y0: vertical coordinate of a point on the splitting / transition line [m]
    :param _ang: rotation angle of the splitting / transition line [rad]
    :param _xc: horizontal coordinate of center [m] of the interval / range [m]
    :param _yc: vertical coordinate of center [m] of the interval / range [m]
    :param _nx: number of points vs horizontal position to represent the transmission element
    :param _ny: number of points vs vertical position to represent the transmission element
    :return: transmission (SRWLOptT) type optical element simulating transition region cutting transverse plane into two half-planes with constant refraction / absorption properties
    """

    cosAng = 1
    sinAng = 0
    if(_ang != 0):
        cosAng = cos(_ang)
        sinAng = sin(_ang)

    xStep = _rx/(_nx - 1)
    yStep = _ry/(_ny - 1)

    halfRx = 0.5*_rx
    halfRy = 0.5*_ry
    halfW = 0.5*_w
    
    ampTr1 = exp(-0.5*_thick/_atten_len1)
    ampTr2 = exp(-0.5*_thick/_atten_len2)
    optPathDif1 = -_delta1*_thick
    optPathDif2 = -_delta2*_thick
    
    op = SRWLOptT(_nx, _ny, _rx, _ry, _x=_xc, _y=_yc, _extTr=1, _alloc_base=[1,0])
    arTr = op.arTr

    ofst = 0
    y = _yc - halfRy
    for iy in range(_ny):

        if(_ang == 0):
            ampTr = 1
            optPathDif = 0

            if(y < _y0 - halfW):
                ampTr = ampTr1
                optPathDif = optPathDif1
            elif(y < _y0 + halfW):
                yRel = y - (_y0 - halfW)
                coef = yRel/_w
                if(_tr_typ == 1): #linear transition
                    ampTr = ampTr1 + coef*(ampTr2 - ampTr1)
                    optPathDif = optPathDif1 + coef*(optPathDif2 - optPathDif1)
            else:
                ampTr = ampTr2
                optPathDif = optPathDif2
                
            for ix in range(_nx):
                arTr[ofst] = ampTr
                arTr[ofst + 1] = optPathDif
                ofst += 2
                    
        else: #_ang != 0
            
            x = _xc - halfRx
            for ix in range(_nx):

                x_mi_x0 = x - _x0
                y_mi_y0 = y - _y0
                #xLoc = _x0 + x_mi_x0*cosAng + y_mi_y0*sinAng #check sign before sinAng and treatment of _x0, _y0
                yLoc = _y0 + y_mi_y0*cosAng - x_mi_x0*sinAng

                if(yLoc < _y0 - halfW):
                    ampTr = ampTr1
                    optPathDif = optPathDif1
                elif(yLoc < _y0 + halfW):
                    yRel = yLoc - (_y0 - halfW) #??
                    coef = yRel/_w
                    if(_tr_typ == 1): #linear transition
                        ampTr = ampTr1 + coef*(ampTr2 - ampTr1)
                        optPathDif = optPathDif1 + coef*(optPathDif2 - optPathDif1)
                else:
                    ampTr = ampTr2
                    optPathDif = optPathDif2

                arTr[ofst] = ampTr
                arTr[ofst + 1] = optPathDif
                ofst += 2

                x += xStep
        y += yStep
    return op

#****************************************************************************
def srwl_opt_setup_gen_transm(_func_path, _delta, _atten_len, _rx, _ry, _xc=0, _yc=0, _ext_tr=0, _fx=0, _fy=0, _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
    """
    Setup Transmission type Optical Element similar to the one simulating CRL, but with arbitrary optical path in material over hor. and vert. positions, defined by external function _func_path(x, y)
    :param _func_path: user-defined function of 2 variables specifying path in material over hor. and vert. positions (x, y)
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _rx: horizontal aperture (range) size [m]
    :param _ry: vertical aperture (range) size [m]
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :param _ext_tr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1)
    :param _fx: horizontal focal length [m]; if it is not set, a numerical estimate will be made (around _xc, _yc)
    :param _fy: vertical focal length [m]; if it is not set, a numerical estimate will be made (around _xc, _yc)
    :param _e_start: initial photon energy
    :param _e_fin: final photon energy
    :param _nx: number of points vs horizontal position to represent the transmission element
    :param _ny: number of points vs vertical position to represent the transmission element
    :return: transmission (SRWLOptT) type optical element which simulates CRL
    """

    ne = 1
    arDelta = [0]
    arAttenLen = [0]

    #if(((type(_delta).__name__ == 'list') or (type(_delta).__name__ == 'array')) and ((type(_atten_len).__name__ == 'list') or (type(_atten_len).__name__ == 'array'))):
    if(isinstance(_delta, list) or isinstance(_delta, array)) and (isinstance(_atten_len, list) or isinstance(_atten_len, array)):
        ne = len(_delta)
        ne1 = len(_atten_len)
        if(ne > ne1):
            ne = ne1
        arDelta = _delta
        arAttenLen = _atten_len
    else:
        arDelta[0] = _delta
        arAttenLen[0] = _atten_len

    dx = _rx/(_nx - 1)
    dy = _ry/(_ny - 1)

    nTot = 2*ne*_nx*_ny #total array length to store amplitude transmission and optical path difference
    arLocTr = array('d', [0]*nTot)
    #Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    ofst = 0
    y = -0.5*_ry #Transmission is always centered on the grid, however grid can be shifted (?)
    for iy in range(_ny):
        x = -0.5*_rx
        for ix in range(_nx):
            pathDif = _func_path(x, y)
            for ie in range(ne):
                arLocTr[ofst] = exp(-0.5*pathDif/arAttenLen[ie]) #amplitude transmission
                arLocTr[ofst + 1] = -arDelta[ie]*pathDif #optical path difference
                ofst += 2
            x += dx
        y += dy

    #Estimating focal lengths
    fx = _fx
    fy = _fy
    avgDelta = arDelta[int(0.5*ne)]
    if(avgDelta != 0):
        if(fx == 0):
            fm1 = _func_path(_xc - dx, _yc)
            f0 = _func_path(_xc, _yc)
            f1 = _func_path(_xc + dx, _yc)
            dfdx = 0.5*(f1 - fm1)/dx
            d2fdx2 = (fm1 - 2*f0 + f1)/(dx*dx)
            fx = 1.e+23
            if(d2fdx2 != 0):
                auxSqrt = sqrt(1 + dfdx*dfdx)
                radCurvX = auxSqrt*auxSqrt*auxSqrt/d2fdx2
                #print('Estimated radCurvX:', radCurvX)
                fx = radCurvX/avgDelta #see the formula for CRL
                #print('Estimated Fx=', fx)
        if(fy == 0):
            fm1 = _func_path(_xc, _yc - dy)
            f0 = _func_path(_xc, _yc)
            f1 = _func_path(_xc, _yc + dy)
            dfdy = 0.5*(f1 - fm1)/dy
            d2fdy2 = (fm1 - 2*f0 + f1)/(dy*dy)
            fy = 1.e+23
            if(d2fdy2 != 0):
                auxSqrt = sqrt(1 + dfdy*dfdy)
                radCurvY = auxSqrt*auxSqrt*auxSqrt/d2fdy2
                fy = radCurvY/avgDelta #see the formula for CRL
                #print('Estimated Fy=', fy)

    return SRWLOptT(_nx, _ny, _rx, _ry, arLocTr, _ext_tr, fx, fy, _xc, _yc, ne, _e_start, _e_fin)

#****************************************************************************
def srwl_opt_setup_zp_tilt(_zp, _ang_tilt, _ang_ax=1.5707963267949, _n_slices=1, _delta=None, _atten_len=None, _e_start=0, _e_fin=0, _nx=None, _ny=None):
    """
    Under Development!!!
    Setup Contained of Transmission type Optical Elements and Drift Spaces simulating tilted thick Zone Plate
    :param _zp: Zone Plate object (of type SRWLOptZP)
    :param _ang_tilt: tilt angle of the ZP [rad]
    :param _ang_ax: angle [rad] from horizontal axis defining orientation of ZP rotation / tilt axis (i.e. if it is equal to pi/2 zp is assumed to be rotated about vertical axis in the horizontal plane)
    :param _n_slices: number of "slices" in longitudinal direction to be used to represent the tilted "thick" Zone Plate
    :param _e_start: initial photon energy [eV]
    :param _e_fin: final photon energy [eV]
    :param _nx: number of points vs horizontal position to represent the transmission element
    :param _ny: number of points vs vertical position to represent the transmission element
    :return: optical element container (of type SRWLOptC) with Transmission optical elements (of type SRWLOptT) separated by drift spaces (of type SRWLOptD)
    """

    if((_zp is None) or (not isinstance(_zp, SRWLOptZP))): raise Exception("Zone Plate object was not submitted.")

    #Set _nx, _ny for Transmission objects from ZP params
    nxSetupReq = False
    if(_nx is None): nxSetupReq = True
    if(_nx <= 0): nxSetupReq = True
    nySetupReq = False
    if(_ny is None): nySetupReq = True
    if(_ny <= 0): nySetupReq = True

    nx = 10*_zp.nZones if nxSetupReq else _nx #based on analytical estimate, 4*k*n, k~=3
    ny = 10*_zp.nZones if nySetupReq else _ny #based on analytical estimate, 4*k*n, k~=3

    #ne = 1
    #arDelta = [0]
    #arAttenLen = [0]

    #if(isinstance(_delta, list) or isinstance(_delta, array)) and (isinstance(_atten_len, list) or isinstance(_atten_len, array)):
    #    ne = len(_delta)
    #    ne1 = len(_atten_len)
    #    if(ne > ne1):
    #        ne = ne1
    #    arDelta = _delta
    #    arAttenLen = _atten_len
    #else:
    #    arDelta[0] = _delta
    #    arAttenLen[0] = _atten_len

    #Setup the object in C++ function

    return None

#****************************************************************************
class SRWLDet(object):
    """Detector of Radiation"""

    def __init__(self, _xStart=0, _xFin=0, _nx=1, _yStart=0, _yFin=0, _ny=1, _dx=0, _dy=0, _spec_eff=1, _eStart=0, _eFin=0, _ord_interp=1): #, _ne=1): #, _eff_e=False):
        """
        :param _xStart: initial value of horizontal position [m] (/angle [rad])
        :param _xFin: final value of horizontal position [m] (/angle [rad])
        :param _nx: number of points vs horizontal position [m] (/angle [rad])
        :param _yStart: initial value of vertical position [m] (/angle [rad])
        :param _yFin: final value of vertical position [m] (/angle [rad])
        :param _ny: number of points vs vertical position [m] (/angle [rad])
        :param _dx: horizontal pixel size [m] (/angle [rad])
        :param _dy: vertical pixel size [m] (/angle [rad])
        :param _spec_eff: one number or array defining detector spectral efficiency (as function of photon energy, on the mesh given by _eStart, _eFin, _ne)
        :param _eStart: initial value of photon energy [eV] (/time [s])
        :param _eFin: final value of photon energy [eV] (/time [s])
        :param _ne: number of points vs photon energy [eV] (/time [s])
        :param _eff_e: treat spectral efficiency as for electric field if True
        :param _ord_interp: interpolation order (i.e. order of polynomials to be used at 2D interpolation)
        """
        self.xStart = _xStart
        self.xFin = _xFin
        self.nx = _nx
        self.yStart = _yStart
        self.yFin = _yFin
        self.ny = _ny
        self.dx = _dx
        self.dy = _dy
        self.specEff = _spec_eff
        self.eStartEff = _eStart
        self.eFinEff = _eFin
        self.ord_interp = _ord_interp

        ne = 1
        if(isinstance(_spec_eff, list) or isinstance(_spec_eff, array)): ne = len(_spec_eff)
        self.neEff = ne

        self.eStepEff = 0 if(ne <= 1) else (_eFin - _eStart)/ne

    def treat_int(self, _stk, _mesh=None, _ord_interp=0):
        """Treat Intensity of Input Radiation
        :param _stk: Stokes data structure or a simple Intensity array to be treated by detector
        :param _mesh: mesh structure (if it is defined, _stk should be treated as Intensity array of type f)
        :param _ord_interp: Interpolation order to be used (overrides self.ord_interp)
        """
        #Move this function to C in the future?

        resInt = SRWLStokes(1, 'f', self.eStartEff, self.eFinEff, 1, self.xStart, self.xFin, self.nx, self.yStart, self.yFin, self.ny, _n_comp=1) #OC20082021 (?)
        #resInt = SRWLStokes(1, 'f', self.eStartEff, self.eFinEff, 1, self.xStart, self.xFin, self.nx, self.yStart, self.yFin, self.ny)

        meshIn = None
        sktIn = None
        if(_mesh is not None) and (isinstance(_mesh, SRWLRadMesh)): 
            meshIn = _mesh
            arI = None
            if(isinstance(_stk, array) and (_stk.itemsize == 4)): #array of type f
                arI = _stk
            else:
                nTot = meshIn.ne*meshIn.nx*meshIn.ny
                arI = array('f', [0]*nTot)
                for i in range(nTot): arI[i] = _stk[i]

            #DEBUG
            #print('treat_int: meshIn.xStart=', meshIn.xStart, ' meshIn.xFin=', meshIn.xFin)
            
            sktIn = SRWLStokes(arI, 'f', meshIn.eStart, meshIn.eFin, meshIn.ne, meshIn.xStart, meshIn.xFin, meshIn.nx, meshIn.yStart, meshIn.yFin, meshIn.ny, _n_comp=1) #OC20082021 (?)
            #sktIn = SRWLStokes(arI, 'f', meshIn.eStart, meshIn.eFin, meshIn.ne, meshIn.xStart, meshIn.xFin, meshIn.nx, meshIn.yStart, meshIn.yFin, meshIn.ny)
            
        else: 
            if(not isinstance(_stk, SRWLStokes)): raise Exception('An object of SRWLStokes class is expected')
            meshIn = _stk.mesh
            #extMeshIsDef = True

        eRange = meshIn.eFin - meshIn.eStart
        eStep = 0 if(meshIn.ne <= 1) else eRange/(meshIn.ne - 1)
        bwMult = 1 if(eRange == 0) else 1./eRange
        ePh = meshIn.eStart
        for ie in range(meshIn.ne):
            effMult = 1 #To treat Spectral Efficiency
            if(self.specEff is not None):
                if(isinstance(self.specEff, list) or isinstance(self.specEff, array)):
                    effMult = uti_math.interp_1d(ePh, self.eStartEff, self.eStepEff, self.neEff, self.specEff, _ord=2) #OC04102021
                    #effMult = interp_1d(ePh, self.eStartEff, self.eStepEff, self.neEff, self.specEff, _ord=2)
                else: effMult = self.specEff

            if((self.dx <= 0) or (self.dy <= 0)):
                ordInterp = _ord_interp if(_ord_interp > 0) else self.ord_interp
                resInt.avg_update_interp(sktIn, _iter=0, _ord=ordInterp, _n_stokes_comp=1, _mult=effMult*bwMult) #to treat all Stokes components / Polarization in the future
            #else: 
                #Program integration within pixels self.dx, self.dy here
            ePh += eStep

        return resInt

    def get_mesh(self):
        eAvg = 0.5*(self.eStartEff + self.eFinEff) #?
        return SRWLRadMesh(
            _eStart=eAvg, _eFin=eAvg, _ne=1, 
            _xStart=self.xStart, _xFin=self.xFin, _nx=self.nx, 
            _yStart=self.yStart, _yFin=self.yFin, _ny=self.ny)

#****************************************************************************
#****************************************************************************
#Auxiliary utility functions
#****************************************************************************
#****************************************************************************
#Moved to uti_math.py:
##def srwl_uti_interp_1d(_x, _x_min, _x_step, _nx, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
##def srwl_uti_interp_2d(_x, _y, _x_min, _x_step, _nx, _y_min, _y_step, _ny, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):

#****************************************************************************
def srwl_uti_ph_en_conv(_x, _in_u='keV', _out_u='nm'):
    """Photon Energy <-> Wavelength conversion
    :param _x: value to be converted
    :param _in_u: input unit
    :param _out_u: output unit
    :return: value in the output units
    """
    #convert _in_u -> [keV]:
    x_keV = _x
    if _in_u == 'eV': x_keV *= 1.e-03
    elif _in_u == '1/cm': x_keV *= (_Light_eV_mu*1.e-07)
    elif _in_u == 'A': x_keV = 10*_Light_eV_mu/x_keV
    elif _in_u == 'nm': x_keV = _Light_eV_mu/x_keV
    elif _in_u == 'um': x_keV = (_Light_eV_mu*1.e-03)/x_keV #this had to be modofoed because of non-ascii "mu" symbol that did not compile on Py 2.7
    elif _in_u == 'mm': x_keV = (_Light_eV_mu*1.e-06)/x_keV
    elif _in_u == 'm': x_keV = (_Light_eV_mu*1.e-09)/x_keV
    elif _in_u == 'THz': x_keV *= (_Light_eV_mu*1000./_LightSp) #sinp="THz";outputval=val*(4.1356672e-06)
    #convert [keV] -> _out_u:
    x = x_keV
    if _out_u == 'eV': x *= 1000.
    elif _out_u == '1/cm': x /= (_Light_eV_mu*1.e-07)
    elif _out_u == 'A': x = 10*_Light_eV_mu/x
    elif _out_u == 'nm': x = _Light_eV_mu/x
    elif _out_u == 'um': x = (_Light_eV_mu*1.e-03)/x #this had to be modifoed because of non-ascii "mu" symbol that did not compile on Py 2.7
    elif _out_u == 'mm': x = (_Light_eV_mu*1.e-06)/x
    elif _out_u == 'm': x = (_Light_eV_mu*1.e-09)/x
    elif _out_u == 'THz': x /= (_Light_eV_mu*1000./_LightSp) #sout="THz";outputval=outputval/(4.1356672e-06)
    return x

#****************************************************************************
def srwl_uti_num_round(_x, _ndig=8):
    order = round(log10(_x))
    fact = 10**order
    return round(_x/fact, _ndig)*fact

#****************************************************************************
def srwl_uti_rand_fill_vol(_np, _x_min, _x_max, _nx, _ar_y_vs_x_min, _ar_y_vs_x_max, _y_min, _y_max, _ny, _ar_z_vs_xy_min, _ar_z_vs_xy_max):
    """
    Generate coordinates of ponts randomly filling 3D volume limited by two arbitrary curves (defining base) and two surfaces
    :param _np: number of random points in rectangular parallelepiped to try
    :param _x_min: min. x coordinate
    :param _x_max: max. x coordinate
    :param _nx: number of points vs x coord.
    :param _ar_y_vs_x_min: min. y vs x array
    :param _ar_y_vs_x_max: max. y vs x array
    :param _y_min: min. y coordinate
    :param _y_max: max. y coordinate
    :param _ny: number of points vs y coord.
    :param _ar_z_vs_xy_min: min. z vs x and y flat 2D array
    :param _ar_z_vs_xy_max: max. z vs x and y flat 2D array
    :return: flat array of point coordinates: array('d', [x1,y1,z1,x2,y2,z2,...])
    """
    yMin = _ar_y_vs_x_min[0]
    yMax = _ar_y_vs_x_max[0]
    for ix in range(_nx):
        yMinCur = _ar_y_vs_x_min[ix]
        if(yMin > yMinCur):
            yMin = yMinCur
        yMaxCur = _ar_y_vs_x_max[ix]
        if(yMax < yMaxCur):
            yMax = yMaxCur

    nxy = _nx*_ny
    zMin = _ar_z_vs_xy_min[0]
    zMax = _ar_z_vs_xy_max[0]
    for ixy in range(nxy):
        zMinCur = _ar_z_vs_xy_min[ixy]
        if(zMin > zMinCur):
            zMin = zMinCur
        zMaxCur = _ar_z_vs_xy_max[ixy]
        if(zMax < zMaxCur):
            zMax = zMaxCur

    xStep = (_x_max - _x_min)/(_nx - 1)
    yStep = (_y_max - _y_min)/(_ny - 1)
    xCen = 0.5*(_x_min + _x_max)
    yCen = 0.5*(yMin + yMax)
    zCen = 0.5*(zMin + zMax)
    xRange = _x_max - _x_min
    yRange = yMax - yMin
    zRange = zMax - zMin

    arPtCoord = array('d', [0]*(_np*3))
    iPtCount = 0
    random.seed()
    for i in range(_np):
        x = xCen + xRange*(random.random() - 0.5)
        y = yCen + yRange*(random.random() - 0.5)
        #yTestMin = srwl_uti_interp_1d(x, _x_min, xStep, _nx, _ar_y_vs_x_min)
        yTestMin = uti_math.interp_1d(x, _x_min, xStep, _nx, _ar_y_vs_x_min)
        #yTestMax = srwl_uti_interp_1d(x, _x_min, xStep, _nx, _ar_y_vs_x_max)
        yTestMax = uti_math.interp_1d(x, _x_min, xStep, _nx, _ar_y_vs_x_max)
        
        if((y >= yTestMin) and (y <= yTestMax)):
            z = zCen + zRange*(random.random() - 0.5)
            #zTestMin = srwl_uti_interp_2d(x, y, _x_min, xStep, _nx, _y_min, yStep, _ny, _ar_z_vs_xy_min)
            zTestMin = uti_math.interp_2d(x, y, _x_min, xStep, _nx, _y_min, yStep, _ny, _ar_z_vs_xy_min)
            #zTestMax = srwl_uti_interp_2d(x, y, _x_min, xStep, _nx, _y_min, yStep, _ny, _ar_z_vs_xy_max)
            zTestMax = uti_math.interp_2d(x, y, _x_min, xStep, _nx, _y_min, yStep, _ny, _ar_z_vs_xy_max)
            if((z >= zTestMin) and (z <= zTestMax)):
                ofst = iPtCount*3
                arPtCoord[ofst] = x
                arPtCoord[ofst + 1] = y
                arPtCoord[ofst + 2] = z
                iPtCount += 1

    if(iPtCount == _np):
        return arPtCoord
    else: #is there faster way to truncate array?
        nResCoord = iPtCount*3
        arResPtCoord = array('d', [0]*nResCoord)
        for i in range(nResCoord):
            arResPtCoord[i] = arPtCoord[i]
            
        return arResPtCoord

#****************************************************************************
def srwl_uti_proc_is_master():
    """
    Check if process is Master (in parallel processing sense)
    """
    try:
        ##resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], -1) #MPI module dynamic load
        #resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], 0) #MPI module dynamic load
        ##multiple re-import won't hurt; but it would be better to avoid this(?)
        #MPI = resImpMPI4Py.MPI
        
        from mpi4py import MPI #OC091014
        
        comMPI = MPI.COMM_WORLD
        rankMPI = comMPI.Get_rank()
        if(rankMPI == 0):
            return True
        else:
            return False
    except:
        return True

#**********************Auxiliary function to write tabulated resulting Intensity data to an ASCII file:
def srwl_uti_save_intens_ascii(_ar_intens, _mesh, _file_path, _n_stokes=1, _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'], _mutual=0, _cmplx=0): #OC06052018
#def srwl_uti_save_intens_ascii(_ar_intens, _mesh, _file_path, _n_stokes=1, _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'], _mutual=0):
    f = open(_file_path, 'w')
    arLabelUnit = [_arLabels[i] + ' [' + _arUnits[i] + ']' for i in range(4)]

    sUnitEnt = arLabelUnit[3]

    #if(_mutual != 0):
    if(_mutual == 1): #OC16072019 (to allow e.g. _mutual == 2 when header should not be modified)
        sUnitEntParts = ['']
        if((sUnitEnt is not None) and (len(sUnitEnt) > 0)): sUnitEntParts = sUnitEnt.split(' ')
        sUnitEntTest = sUnitEnt
        if(len(sUnitEntParts) > 0): sUnitEntTest = sUnitEntParts[0]
        sUnitEntTest = sUnitEntTest.replace(' ', '')
        if(sUnitEntTest.lower != 'mutual'):
            sPrefix = 'Mutual' #this prefix is a switch meaning eventual special processing in viewing utilities
            if(_cmplx != 0): sPrefix = 'Complex Mutual' #OC06052018
            if(sUnitEnt.startswith(' ') == False): sPrefix += ' '
            sUnitEnt = sPrefix + sUnitEnt

    #print(sUnitEnt) #DEBUG
    
    f.write('#' + sUnitEnt + ' (C-aligned, inner loop is vs ' + _arLabels[0] + ', outer loop vs ' + _arLabels[2] + ')\n')
    f.write('#' + repr(_mesh.eStart) + ' #Initial ' + arLabelUnit[0] + '\n')
    f.write('#' + repr(_mesh.eFin) + ' #Final ' + arLabelUnit[0] + '\n')
    f.write('#' + repr(_mesh.ne) + ' #Number of points vs ' + _arLabels[0] + '\n')
    f.write('#' + repr(_mesh.xStart) + ' #Initial ' + arLabelUnit[1] + '\n')
    f.write('#' + repr(_mesh.xFin) + ' #Final ' + arLabelUnit[1] + '\n')
    f.write('#' + repr(_mesh.nx) + ' #Number of points vs ' + _arLabels[1] + '\n')
    f.write('#' + repr(_mesh.yStart) + ' #Initial ' + arLabelUnit[2] + '\n')
    f.write('#' + repr(_mesh.yFin) + ' #Final ' + arLabelUnit[2] + '\n')
    f.write('#' + repr(_mesh.ny) + ' #Number of points vs ' + _arLabels[2] + '\n')

    #strOut =  '#' + sUnitEnt + ' (C-aligned, inner loop is vs ' + _arLabels[0] + ', outer loop vs ' + _arLabels[2] + ')\n'
    #strOut += '#' + repr(_mesh.eStart) + ' #Initial ' + arLabelUnit[0] + '\n'
    #strOut += '#' + repr(_mesh.eFin) + ' #Final ' + arLabelUnit[0] + '\n'
    #strOut += '#' + repr(_mesh.ne) + ' #Number of points vs ' + _arLabels[0] + '\n'
    #strOut += '#' + repr(_mesh.xStart) + ' #Initial ' + arLabelUnit[1] + '\n'
    #strOut += '#' + repr(_mesh.xFin) + ' #Final ' + arLabelUnit[1] + '\n'
    #strOut += '#' + repr(_mesh.nx) + ' #Number of points vs ' + _arLabels[1] + '\n'
    #strOut += '#' + repr(_mesh.yStart) + ' #Initial ' + arLabelUnit[2] + '\n'
    #strOut += '#' + repr(_mesh.yFin) + ' #Final ' + arLabelUnit[2] + '\n'
    #strOut += '#' + repr(_mesh.ny) + ' #Number of points vs ' + _arLabels[2] + '\n'
            
    nComp = 1
    if _n_stokes > 0:
        f.write('#' + repr(_n_stokes) + ' #Number of components\n')
        #DEBUG
        #print('#' + repr(_n_stokes) + ' #Number of components\n')
        
        #strOut += '#' + repr(_n_stokes) + ' #Number of components\n'
        nComp = _n_stokes
    nRadPt = _mesh.ne*_mesh.nx*_mesh.ny
    if(_mutual > 0): nRadPt *= nRadPt
    
    nVal = nRadPt*nComp #_mesh.ne*_mesh.nx*_mesh.ny*nComp
    if(_cmplx != 0): nVal *= 2 #OC06052018

    #DEBUG
    #print('Saving (header) Hor. Mesh:', _mesh.xStart, _mesh.xFin, _mesh.nx)
    #print('Saving (header) Vert. Mesh:', _mesh.yStart, _mesh.yFin, _mesh.ny)
    #print('Saving (header) _mutual:', _mutual)
    #print('Saving (header) nComp:', nComp)
    #END DEBUG

    for i in range(nVal): #write all data into one column using "C-alignment" as a "flat" 1D array
        f.write(' ' + repr(_ar_intens[i]) + '\n')
        #strOut += ' ' + repr(_ar_intens[i]) + '\n'
        #DEBUG
        #if(_ar_intens[i] != 0.): print('i=', i, ' Non-zero Int. value:', _ar_intens[i])
        #if(i > 450000): break
        #END DEBUG
    
    #f = open(_file_path, 'w')
    #f.write(strOut)
    f.close()

#**********************Auxiliary function to read-in tabulated  Intensity data from an ASCII file (format is defined in srwl_uti_save_intens_ascii)
def srwl_uti_read_intens_ascii(_file_path, _num_type='f'):

    sCom = '#'
    f = open(_file_path, 'r')
    lines = f.readlines()

    resMesh = SRWLRadMesh()

    curParts = lines[1].split(sCom); resMesh.eStart = float(curParts[1]) #to check
    curParts = lines[2].split(sCom); resMesh.eFin = float(curParts[1]) #to check
    curParts = lines[3].split(sCom); resMesh.ne = int(curParts[1]) #to check
    
    curParts = lines[4].split(sCom); resMesh.xStart = float(curParts[1]) #to check
    curParts = lines[5].split(sCom); resMesh.xFin = float(curParts[1]) #to check
    curParts = lines[6].split(sCom); resMesh.nx = int(curParts[1]) #to check

    curParts = lines[7].split(sCom); resMesh.yStart = float(curParts[1]) #to check
    curParts = lines[8].split(sCom); resMesh.yFin = float(curParts[1]) #to check
    curParts = lines[9].split(sCom); resMesh.ny = int(curParts[1]) #to check

    iStart = 10
    if((lines[10])[0] == sCom): iStart = 11
    
    nRows = len(lines)
    arInt = []
    for i in range(iStart, nRows):
        curLine = lines[i]
        if(len(curLine) > 0): arInt.append(float(curLine))
    f.close()
    return array(_num_type, arInt), resMesh

#**********************Auxiliary function to write tabulated resulting Intensity data to an HDF5 file:
def srwl_uti_save_intens_hdf5(_ar_intens, _mesh, _file_path, _n_stokes=1,
                              _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'],
                              _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'], _mutual=0, _cmplx=0, _wfr=None): #OC05052021
                              #_arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'], _mutual=0, _cmplx=0): #RAC30032020
    #To review!
    ### Load package Numpy
    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')
        #print('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following ' + 
        #      "command to install it: \npip install numpy")
        
    ### Load package h5py
    try:
        import h5py as h5
    except:
        raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')
        #print('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following ' + 
        #      "command to install it: \npip install h5py")

    ### Begin time record of creating and placing objects
    #t0 = time.time();
    
    #Get argument Lables
    #OC14022021 (commented-out)
    #arLabelUnit = [_arLabels[i] + ' [' + _arUnits[i] + ']' for i in range(4)]
    #arLabelUnit_ascii = [xx.encode("ascii", "ignore") for xx in arLabelUnit] #convert from U19 to ascii for h5py
    #arLabels_ascii = [xx.encode("ascii", "ignore") for xx in _arLabels] #convert from U19 to ascii for h5py
    #sUnitEnt = arLabelUnit[3]

    #OC14022021 (commented-out)
    #if(_mutual != 0):
    #    sUnitEntParts = ['']
    #    if((sUnitEnt is not None) and (len(sUnitEnt) > 0)): sUnitEntParts = sUnitEnt.split(' ')
    #    sUnitEntTest = sUnitEnt
    #    if(len(sUnitEntParts) > 0): sUnitEntTest = sUnitEntParts[0]
    #    sUnitEntTest = sUnitEntTest.replace(' ', '')
    #    if(sUnitEntTest.lower != 'mutual'):
    #        sPrefix = 'Mutual' #this prefix is a switch meaning eventual special processing in viewing utilities
    #        if(_cmplx != 0): sPrefix = 'Complex Mutual'
    #        if(sUnitEnt.startswith(' ') == False): sPrefix += ' '
    #        sUnitEnt = sPrefix + sUnitEnt

    #Create intensity data set as numpy array
    nComp = 1
    if _n_stokes > 0:
        nComp = _n_stokes
    nRadPt = _mesh.ne*_mesh.nx*_mesh.ny
    if(_mutual > 0): nRadPt *= nRadPt
    nVal = nRadPt*nComp #_mesh.ne*_mesh.nx*_mesh.ny*nComp
    if(_cmplx != 0): nVal *= 2

    if(not isinstance(_ar_intens, (array, np.ndarray))): #OC04102021
    #if(not (isinstance(_ar_intens, array) or isinstance(_ar_intens, np.array))):
        intensity_data = np.array([0]*nVal, 'f')
        intensity_data[:] = _ar_intens
        #intensity_data = np.array([_ar_intens[ii] for ii in range(nVal)])
    else: intensity_data = _ar_intens

    #Write file header and data set (compound datatype) as hdf5
    with h5.File(_file_path, 'w') as hf:
        #intensity
        hf.create_dataset('intensity', data=intensity_data) #Change this name?

        #OC14022021
        ats = hf.attrs
        ats.create('eStart', _mesh.eStart)
        ats.create('eFin', _mesh.eFin)
        ats.create('ne', _mesh.ne)
        ats.create('xStart', _mesh.xStart)
        ats.create('xFin', _mesh.xFin)
        ats.create('nx', _mesh.nx)
        ats.create('yStart', _mesh.yStart)
        ats.create('yFin', _mesh.yFin)
        ats.create('ny', _mesh.ny)
        ats.create('n_stokes', _n_stokes)
        ats.create('mutual', _mutual)
        ats.create('cmplx', _cmplx)
        ats.create('arLabels', _arLabels)
        ats.create('arUnits', _arUnits)

        if(hasattr(_mesh, 'itStart')): 
            #print('_mesh.itStart=', _mesh.itStart) #DEBUG
            ats.create('itStart', _mesh.itStart) #OC17042021
        if(hasattr(_mesh, 'itFin')): 
            #print('_mesh.itFin', _mesh.itFin) #DEBUG
            ats.create('itFin', _mesh.itFin) #OC17042021

        if(hasattr(_mesh, 'type')): ats.create('type', _mesh.type) #OC20042021

        #OC05052021
        if(_wfr is not None): #Average wavefront params (to be used e.g. for subsequent CMD on the MI/CSD data stored in 'intensity')
            ats.create('numTypeElFld', _wfr.numTypeElFld)
            ats.create('Rx', _wfr.Rx)
            ats.create('Ry', _wfr.Ry)
            ats.create('dRx', _wfr.dRx)
            ats.create('dRy', _wfr.dRy)
            ats.create('xc', _wfr.xc)
            ats.create('yc', _wfr.yc)
            ats.create('avgPhotEn', _wfr.avgPhotEn)
            ats.create('presCA', _wfr.presCA)
            ats.create('presFT', _wfr.presFT)
            ats.create('unitElFld', _wfr.unitElFld)
            ats.create('unitElFldAng', _wfr.unitElFldAng)
    
    #DEBUG
    #print('Saving HDF5 with _n_stokes=', _n_stokes, ' _mutual=', _mutual, ' _cmplx=', _cmplx)
    #END DEBUG
    
    hf.close()
    #Print competed time
    #print('HDF5 completed (lasted', round(time.time() - t0, 6), 's)')

#**********************Auxiliary function to read-in intensity data from an hdf5 file (format is defined in srwl_uti_save_intens_hdf5)
def srwl_uti_read_intens_hdf5(_file_path, _num_type='f'): #RAC30032020
    #To review!
    ### Load package Numpy
    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')
        #print('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following ' + 
        #      "command to install it: \npip install numpy")
    ### Load package h5py
    try:
        import h5py as h5
    except:
        raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')
        #print('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following ' + 
        #      "command to install it: \npip install h5py")

    hf = h5.File(_file_path, 'r')
    mesh = SRWLRadMesh()

    #Get attributes
    #OC15022021
    ats = hf.attrs
    mesh.eStart = float(ats.get('eStart'))
    mesh.eFin = float(ats.get('eFin'))
    mesh.ne = int(ats.get('ne'))
    mesh.xStart = float(ats.get('xStart'))  
    mesh.xFin = float(ats.get('xFin'))
    mesh.nx = int(ats.get('nx'))
    mesh.yStart = float(ats.get('yStart'))
    mesh.yFin = float(ats.get('yFin'))
    mesh.ny = int(ats.get('ny'))

    n_stokes = int(ats.get('n_stokes'))
    mutual = int(ats.get('mutual'))
    cmplx = int(ats.get('cmplx'))
    arLabels = ats.get('arLabels')
    arUnits = ats.get('arUnits')

    #OC17042021
    itStart = None
    try: itStart = ats.get('itStart')
    except: pass
    if(itStart is not None): mesh.itStart = int(itStart)
    itFin = None
    try: itFin = ats.get('itFin')
    except: pass
    if(itFin is not None): mesh.itFin = int(itFin)

    if(mutual != 0): mesh.type = 2 #OC20042021 (this indicates MI in C++)
    else: mesh.type = 0

    #OC06052021
    numTypeElFld = None
    try: numTypeElFld = ats.get('numTypeElFld')
    except: pass
    
    wfr = None #auxiliary wavefront
    if(numTypeElFld is not None):
        wfr = SRWLWfr()
        wfr.numTypeElFld = numTypeElFld
        wfr.Rx = float(ats.get('Rx'))
        wfr.Ry = float(ats.get('Ry'))
        wfr.dRx = float(ats.get('dRx'))
        wfr.dRy = float(ats.get('dRy'))
        wfr.xc = float(ats.get('xc'))
        wfr.yc = float(ats.get('yc'))
        wfr.avgPhotEn = float(ats.get('avgPhotEn'))
        wfr.presCA = int(ats.get('presCA'))
        wfr.presFT = int(ats.get('presFT'))
        wfr.unitElFld = int(ats.get('unitElFld'))
        wfr.unitElFldAng = int(ats.get('unitElFldAng'))

        wfr.mesh = mesh #OC19052021

    #Get params from hdf5 data sets
    #OC15022021 (commented-out)
    #resMesh.eStart = float(hf.get('eStart'))
    #resMesh.eFin = float(hf.get('eFin'))
    #resMesh.ne = int(hf.get('ne'))
    #resMesh.xStart = float(hf.get('xStart'))  
    #resMesh.xFin = float(hf.get('xFin'))
    #resMesh.nx = int(hf.get('nx'))
    #resMesh.yStart = float(hf.get('yStart'))
    #resMesh.yFin = float(hf.get('yFin'))
    #resMesh.ny = int(hf.get('ny'))

    #Get intensity data from hdf5 data set
    arInt = np.array(hf.get('intensity'), dtype=_num_type) #OC15022021
    #arInt = np.array(hf.get('intensity_data'), dtype=_num_type)

    #DEBUG
    np.asarray_chkfinite(arInt, dtype=_num_type)
    #END DEBUG

    return arInt, mesh, {'n_stokes': n_stokes, 'mutual': mutual, 'cmplx': cmplx, 'arLabels': arLabels, 'arUnits': arUnits}, wfr #OC06052021
    #return arInt, mesh, {'n_stokes': n_stokes, 'mutual': mutual, 'cmplx': cmplx, 'arLabels': arLabels, 'arUnits': arUnits} #OC15022021

#**********************Auxiliary function to write Intensity to a file in a given format (ASCII and HDF5 are currently supported)
def srwl_uti_save_intens(_ar_intens, _mesh, _file_path, _n_stokes=1, 
                         _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'], 
                         _mutual=0, _cmplx=0, _form='ascii', _wfr=None): #OC06052021
                         #_mutual=0, _cmplx=0, _form='ascii'): #OC18022021
    form = _form.upper()
    if((form == 'ASCII') or (form == 'ASC') or (form == 'TXT')): srwl_uti_save_intens_ascii(_ar_intens, _mesh, _file_path, _n_stokes, _arLabels, _arUnits, _mutual, _cmplx) #Add wfr? #OC19072021
    #if((form == 'ASCII') or (form == 'TXT')): srwl_uti_save_intens_ascii(_ar_intens, _mesh, _file_path, _n_stokes, _arLabels, _arUnits, _mutual, _cmplx) #Add wfr?
    elif((form == 'HDF5') or (form == 'H5')): srwl_uti_save_intens_hdf5(_ar_intens, _mesh, _file_path, _n_stokes, _arLabels, _arUnits, _mutual, _cmplx, _wfr) #OC06052021
    #elif((form == 'HDF5') or (form == 'H5')): srwl_uti_save_intens_hdf5(_ar_intens, _mesh, _file_path, _n_stokes, _arLabels, _arUnits, _mutual, _cmplx)

#**********************Auxiliary function to read Intensity from a file in a given format (ASCII and HDF5 are currently supported)
def srwl_uti_read_intens(_file_path, _num_type='f', _form='ascii'): #OC04032021
    form = _form.upper()
    if((form == 'ASCII') or (form == 'ASC') or (form == 'TXT')): return srwl_uti_read_intens_ascii(_file_path, _num_type) #Make these two return exactly same thing? #OC19072021
    #if((form == 'ASCII') or (form == 'TXT')): return srwl_uti_read_intens_ascii(_file_path, _num_type) #Make these two return exactly same thing?
    elif((form == 'HDF5') or (form == 'H5')): return srwl_uti_read_intens_hdf5(_file_path, _num_type)
    else: return None

#**********************Auxiliary function to write Intensity data, simulating Detector data, to an HDF5 file:
def srwl_uti_save_intens_hdf5_exp(_ar_intens, _mesh, _file_path, _exp_type='XPCS', _dt=0., _dist_smp=0., _bm_size_x=10.e-06, _bm_size_y=10.e-06): #OC20082021 (based on format suggestion by M. Rakitin)
    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')
    try:
        import h5py as h5
    except:
        raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')

    if(_exp_type.upper() == 'XPCS'):
        parameters = {
            "time_between_frames": _dt,
            "energy": _mesh.eStart,
            "detector_coordinates": (0.5*(_mesh.xStart + _mesh.xFin), 0.5*(_mesh.yStart + _mesh.yFin)), #Center Position vs X and Y
            "detector_distance_from_sample": _dist_smp,
            "pixel_size_x": (_mesh.xFin - _mesh.xStart)/_mesh.nx, #This is more like distance between neighboring pixel centers
            "pixel_size_y": (_mesh.yFin - _mesh.yStart)/_mesh.ny, #This is more like distance between neighboring pixel centers
            "beam_size_x": _bm_size_x, #Not clear what it is - X-ray beam size at sample?
            "beam_size_y": _bm_size_y,
        }

        with h5.File(_file_path, 'w') as hf:

            # intensities = np.random.random((10, 640, 480))
            group = hf.create_group("srw")
            group.create_dataset("intensities", data=_ar_intens) #_ar_intens should be NumPy array (3D)
            # For maximum compression using gzip, comment out the previous line and uncomment the next one:
            # group.create_dataset('intensities', data=intensities, compression="gzip", compression_opts=9)

            param_subgroup = group.create_group("parameters")
            for k, v in parameters.items():
                param_subgroup.create_dataset(k, data=v)

#**********************Auxiliary function to convert an hdf5 file to an ASCII file
# def srwl_uti_convert_intens_hdf5_to_ascii(_file_path): #RAC30032020
#     #To be re-written!
#     ### Load package Numpy
#     try:
#         import numpy as np
#     except:
#         raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')
#         #print('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following ' + 
#         #      "command to install it: \npip install numpy")
        
#     ### Load package h5py
#     try:
#         import h5py as h5
#     except:
#         raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')
#         #print('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following ' + 
#         #      "command to install it: \npip install h5py")
    
#     hf = h5.File(_file_path, 'r')
    
#     #Get params from hdf5 data sets
#     eStart = float(hf.get('eStart'))
#     eFin = float(hf.get('eFin'))
#     ne = int(hf.get('ne'))
#     xStart = float(hf.get('xStart'))  
#     xFin = float(hf.get('xFin'))
#     nx = int(hf.get('nx'))
#     yStart = float(hf.get('yStart'))
#     yFin = float(hf.get('yFin'))
#     ny = int(hf.get('ny'))
#     n_stokes = int(hf.get('n_stokes'))
#     mutual = int(hf.get('mutual'))
#     cmplx = int(hf.get('cmplx'))
    
#     #Get labels
#     arLabels = hf.get('arLabels')
#     arLabelUnit = hf.get('arLabelUnit')
#     sUnitEnt = hf.get('sUnitEnt')
    
#     #Convert intensity data from hdf5 data set to numpy array
#     intensity_data = hf.get('intensity_data').value #intensity_data is now an ndarray.
    
#     #Create ascii file
#     f = open(_file_path + '-ascii.dat', 'w')
    
#     f.write('#' + sUnitEnt + ' (C-aligned, inner loop is vs ' + arLabels[0] + ', outer loop vs ' + arLabels[2] + ')\n')
#     f.write('#' + repr(eStart) + ' #Initial ' + arLabelUnit[0] + '\n')
#     f.write('#' + repr(eFin) + ' #Final ' + arLabelUnit[0] + '\n')
#     f.write('#' + repr(ne) + ' #Number of points vs ' + arLabels[0] + '\n')
#     f.write('#' + repr(xStart) + ' #Initial ' + arLabelUnit[1] + '\n')
#     f.write('#' + repr(xFin) + ' #Final ' + arLabelUnit[1] + '\n')
#     f.write('#' + repr(nx) + ' #Number of points vs ' + arLabels[1] + '\n')
#     f.write('#' + repr(yStart) + ' #Initial ' + arLabelUnit[2] + '\n')
#     f.write('#' + repr(yFin) + ' #Final ' + arLabelUnit[2] + '\n')
#     f.write('#' + repr(ny) + ' #Number of points vs ' + arLabels[2] + '\n')
            
#     nComp = 1
#     if n_stokes > 0:
#         f.write('#' + repr(n_stokes) + ' #Number of components\n')
#         nComp = n_stokes
#     nRadPt = ne*nx*ny
#     if(mutual > 0): nRadPt *= nRadPt
    
#     nVal = nRadPt*nComp #ne*nx*ny*nComp
#     if(cmplx != 0): nVal *= 2 #OC06052018

#     for ii in range(nVal): #write all data into one column using "C-alignment" as a "flat" 1D array
#         f.write(' ' + repr(intensity_data[ii]) + '\n')

#     f.close()

#**********************Auxiliary function to save Updated Coherent Modes file
def srwl_uti_save_wfr_cm_hdf5(_arEx, _arEy, _awfr, _file_path): #OC06052021
#def srwl_uti_save_wfr_cm_hdf5(_arE, _awfr, _file_path): #OC27062021
#def srwl_uti_save_wfr_cm_hdf5(_file_path, _lst_wfr): #OC06052021
#def srwl_uti_save_wfr_cm_hdf5(_file_path, arE, mesh, _gen0s=True): #RL23022021
    """
    Save CMD data (with a number of fully-coherent wavefronts, calculated in the same mesh vs x nad y) to a file
    :param _arE: 2D array of complex electric field data of the coherent modes
    :param _awfr: wavefront () of mesh 
    :param _file_path: string specifying path do data file to be loaded
    """
    ### Load package numpy
    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')
    ### Load package h5py
    try:
        import h5py as h5
    except:
        raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')

    mesh = None
    wfrDefined = False
    if(isinstance(_awfr, SRWLWfr)): 
        mesh = _awfr.mesh #OC06052021
        wfrDefined = True
    elif(isinstance(_awfr, SRWLRadMesh)): mesh = _awfr #OC18052021
    else: raise Exception('Mesh data was not supplied.')

    #mesh = _wfr0.mesh #OC06052021

    #hf = h5.File(_file_path, 'w')
    with h5.File(_file_path, 'w') as hf:

        #if(_arE is not None): #OC27062021
        if(_arEx is not None): 
            #data = _arE.view(np.float32)
            #hf.create_dataset('arE', data=data, dtype=np.float32, shape=data.shape) #OC27062021
            data = _arEx.view(np.float32)
            hf.create_dataset('arEx', data=data, dtype=np.float32, shape=data.shape) #OC08052021
        if(_arEy is not None): 
            data = _arEy.view(np.float32)
            hf.create_dataset('arEy', data=data, dtype=np.float32, shape=data.shape)
        #if(_arEx is not None): hf.create_dataset('arEx', data=_arEx)
        #if(_arEy is not None): hf.create_dataset('arEy', data=_arEy)

        #Write attributes
        ats = hf.attrs
        ats.create('eStart', mesh.eStart)
        ats.create('eFin', mesh.eFin)
        ats.create('ne', mesh.ne)
        ats.create('xStart', mesh.xStart)
        ats.create('xFin', mesh.xFin)
        ats.create('nx', mesh.nx)
        ats.create('yStart', mesh.yStart)
        ats.create('yFin', mesh.yFin)
        ats.create('ny', mesh.ny)
        ats.create('zStart', mesh.zStart)

        if(wfrDefined):
            ats.create('numTypeElFld', _awfr.numTypeElFld)
            ats.create('Rx', _awfr.Rx)
            ats.create('Ry', _awfr.Ry)
            ats.create('dRx', _awfr.dRx)
            ats.create('dRy', _awfr.dRy)
            ats.create('xc', _awfr.xc)
            ats.create('yc', _awfr.yc)
            ats.create('avgPhotEn', _awfr.avgPhotEn)
            ats.create('presCA', _awfr.presCA)
            ats.create('presFT', _awfr.presFT)
            ats.create('unitElFld', _awfr.unitElFld)
            ats.create('unitElFldAng', _awfr.unitElFldAng)
        else:
            ats.create('numTypeElFld', 'f')
            ats.create('Rx', 0.)
            ats.create('Ry', 0.)
            ats.create('dRx', 0.)
            ats.create('dRy', 0.)
            ats.create('xc', 0.)
            ats.create('yc', 0.)
            ats.create('avgPhotEn', 0.5*(mesh.eStart + mesh.xFin))
            ats.create('presCA', 0)
            ats.create('presFT', 0)
            ats.create('unitElFld', 0) #or 1 - sqrt(Phot/s/0.1%bw/mm^2) ?
            ats.create('unitElFldAng', 0)
            
    #hf.attrs['eStart'] = mesh.eStart #= float(ats.get('eStart'))
    #hf.attrs['eFin'] = mesh.eFin #= float(ats.get('eFin'))
    #hf.attrs['ne'] = mesh.ne# = int(ats.get('ne'))
    #hf.attrs['xStart'] = mesh.xStart #= float(ats.get('xStart'))
    #hf.attrs['xFin'] = mesh.xFin #= float(ats.get('xFin'))
    #hf.attrs['nx'] = mesh.nx #= int(ats.get('nx'))
    #hf.attrs['yStart'] = mesh.yStart# = float(ats.get('yStart'))
    #hf.attrs['yFin'] = mesh.yFin #= float(ats.get('yFin'))
    #hf.attrs['ny'] = mesh.ny #= int(ats.get('ny'))
    #hf.attrs['zStart'] = mesh.zStart #= float(ats.get('zStart'))

    # wfr.numTypeElFld = ats.get('numTypeElFld')
    # wfr.Rx = float(ats.get('Rx'))
    # wfr.Ry = float(ats.get('Ry'))
    # wfr.dRx = float(ats.get('dRx'))
    # wfr.dRy = float(ats.get('dRy'))
    # wfr.xc = float(ats.get('xc'))
    # wfr.yc = float(ats.get('yc'))
    # wfr.avgPhotEn = float(ats.get('avgPhotEn'))
    # wfr.presCA = int(ats.get('presCA'))
    # wfr.presFT = int(ats.get('presFT'))
    # wfr.unitElFld = int(ats.get('unitElFld'))
    # wfr.unitElFldAng = int(ats.get('unitElFldAng'))

    # if False:
    #     wfr.numTypeElFld = ats.get('numTypeElFld')
    #     wfr.Rx = float(ats.get('Rx'))
    #     wfr.Ry = float(ats.get('Ry'))
    #     wfr.dRx = float(ats.get('dRx'))
    #     wfr.dRy = float(ats.get('dRy'))
    #     wfr.xc = float(ats.get('xc'))
    #     wfr.yc = float(ats.get('yc'))
    #     wfr.avgPhotEn = float(ats.get('avgPhotEn'))
    #     wfr.presCA = int(ats.get('presCA'))
    #     wfr.presFT = int(ats.get('presFT'))
    #     wfr.unitElFld = int(ats.get('unitElFld'))
    #     wfr.unitElFldAng = int(ats.get('unitElFldAng'))

    #arEid = ['arEx', 'arEy']
    #Save All Electric Field data sets
    #for i in range(len(arE)):
    #    data = arE[i].view(np.float32)
    #    hf.create_dataset(arEid[i], data=data, dtype='f', shape=data.shape)
    #hf.close()

    #Save All Electric Field data sets
    # if arE[0] is not None: # RL05052021
    #     data = arE[0].view(np.float32)
    #     hf.create_dataset("arEx", data=data, dtype=np.float32, shape=data.shape)
    # if arE[1] is not None:
    #     data = arE[1].view(np.float32)
    #     hf.create_dataset("arEy", data=data, dtype=np.float32, shape=data.shape)

    hf.close()
    return 

#**********************Auxiliary function to read-in Updated Coherent Modes file
def srwl_uti_read_wfr_cm_hdf5(_file_path, _gen0s=True): #OC11042020
    """
    Reads-in Wavefront data file (with a number of wavefronts, calculated in the same mesh vs x nad y)
    :param _file_path: string specifying path do data file to be loaded
    :param _gen0s: switch defining whether zero-electric field array(s) need to be created if some polarization component(s) are missing in the file 
    :return: list of wavefronts (objects of SRWLWfr type)
    """
    ### Load package numpy
    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')

    ### Load package h5py
    try:
        import h5py as h5
    except:
        raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')

    wfr = SRWLWfr() #auxiliary wavefront
    mesh = wfr.mesh

    hf = h5.File(_file_path, 'r')

    #Get attributes
    ats = hf.attrs
    mesh.eStart = float(ats.get('eStart'))
    mesh.eFin = float(ats.get('eFin'))
    mesh.ne = int(ats.get('ne'))
    mesh.xStart = float(ats.get('xStart'))  
    mesh.xFin = float(ats.get('xFin'))
    mesh.nx = int(ats.get('nx'))
    mesh.yStart = float(ats.get('yStart'))
    mesh.yFin = float(ats.get('yFin'))
    mesh.ny = int(ats.get('ny'))
    mesh.zStart = float(ats.get('zStart'))

    #OC02082022
    nMom = 11*mesh.ne
    wfr.arMomX = array('d', [0] * nMom)
    wfr.arMomY = array('d', [0] * nMom)

    try: #OC09052021 (to walk around cases when this input is absent)
        wfr.numTypeElFld = ats.get('numTypeElFld')
        wfr.Rx = float(ats.get('Rx'))
        wfr.Ry = float(ats.get('Ry'))
        wfr.dRx = float(ats.get('dRx'))
        wfr.dRy = float(ats.get('dRy'))
        wfr.xc = float(ats.get('xc'))
        wfr.yc = float(ats.get('yc'))
        wfr.avgPhotEn = float(ats.get('avgPhotEn'))
        wfr.presCA = int(ats.get('presCA'))
        wfr.presFT = int(ats.get('presFT'))
        wfr.unitElFld = int(ats.get('unitElFld'))
        wfr.unitElFldAng = int(ats.get('unitElFldAng'))
    except:
        wfr.numTypeElFld = 'f'
        wfr.Rx = 0
        wfr.Ry = 0
        wfr.dRx = 0
        wfr.dRy = 0
        wfr.xc = 0
        wfr.yc = 0
        wfr.avgPhotEn = 0
        wfr.presCA = 0
        wfr.presFT = 0
        wfr.unitElFld = 1
        wfr.unitElFldAng = 0

    #Get All Electric Field data sets
    arEx = None
    arExH5 = hf.get('arEx')
    if(arExH5 is not None): arEx = np.array(arExH5)

    arEy = None
    arEyH5 = hf.get('arEy')
    if(arEyH5 is not None): arEy = np.array(arEyH5)

    nWfr = 0
    lenArE = 0
    if(arEx is not None):
        nWfr = len(arEx)
        lenArE = len(arEx[0])
    elif(arEy is not None):
        nWfr = len(arEy)
        lenArE = len(arEy[0])

    arE0s = None #OC28062021
    if(_gen0s and (lenArE > 0)): arE0s = np.array([0]*lenArE, 'f') #OC28062021
    #arE0s = None if(lenArE <= 0) else np.array([0]*lenArE, 'f')
    
    lstWfr = []
    for i in range(nWfr):
        newWfr = deepcopy(wfr)
        newWfr.mesh = copy(mesh) #OC07112020
        #newWfr.mesh = mesh

        if(arEx is not None):
            newWfr.arEx = copy(arE0s)
            newWfr.arEx[:] = arEx[i]
        elif(_gen0s): newWfr.arEx = arE0s #maybe duplicate it?
        
        if(arEy is not None):
            newWfr.arEy = copy(arE0s)
            newWfr.arEy[:] = arEy[i]
        elif(_gen0s): newWfr.arEy = arE0s #maybe duplicate it?

        lstWfr.append(newWfr)

    #Tests
    #print(len(lstWfr))
    #print(lstWfr[0].arEx)
    #print(lstWfr[0].arEy)

    return lstWfr

#**********************Auxiliary function to save a Wavefront to a file
def srwl_uti_save_wfr_hdf5(_wfr, _file_path, _form=0): #OC09032024
    """
    Save CMD data (with a number of fully-coherent wavefronts, calculated in the same mesh vs x nad y) to a file
    :param _wfr: wavefront (SRWLWfr) type object
    :param _file_path: string specifying path do data file to be loaded
    :param _form: format to store electric field data (=0: standard alignment in SRW, =1: with photon energy being outmost cycle, ...)
    """
    ### Load package numpy
    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')
    ### Load package h5py
    try:
        import h5py as h5
    except:
        raise Exception('h5py can not be loaded. You may need to install h5py. If you are using pip, you can use the following command to install it: \npip install h5py')

    if(not isinstance(_wfr, SRWLWfr)): 
        raise Exception('Wavefront (SRWLWfr) type object is expected')

    mesh = _wfr.mesh
    if(not isinstance(mesh, SRWLRadMesh)): 
        raise Exception('Wavefront mesh (SRWLRadMesh) type object is expected')
    
    ne = mesh.ne; nx = mesh.nx; ny = mesh.ny
    nxny = nx*ny; nenxny = ne*nxny
    arECR = None

    # if(not isinstance(_ar_intens, (array, np.ndarray))): #OC04102021
    # #if(not (isinstance(_ar_intens, array) or isinstance(_ar_intens, np.array))):
    #     intensity_data = np.array([0]*nVal, 'f')
    #     intensity_data[:] = _ar_intens
    #     #intensity_data = np.array([_ar_intens[ii] for ii in range(nVal)])
    # else: intensity_data = _ar_intens

    with h5.File(_file_path, 'w') as hf:

        #Write attributes first
        ats = hf.attrs
        ats.create('eStart', mesh.eStart)
        ats.create('eFin', mesh.eFin)
        ats.create('ne', ne)
        ats.create('xStart', mesh.xStart)
        ats.create('xFin', mesh.xFin)
        ats.create('nx', nx)
        ats.create('yStart', mesh.yStart)
        ats.create('yFin', mesh.yFin)
        ats.create('ny', ny)
        ats.create('zStart', mesh.zStart)

        ats.create('numTypeElFld', _wfr.numTypeElFld)
        ats.create('Rx', _wfr.Rx)
        ats.create('Ry', _wfr.Ry)
        ats.create('dRx', _wfr.dRx)
        ats.create('dRy', _wfr.dRy)
        ats.create('xc', _wfr.xc)
        ats.create('yc', _wfr.yc)
        ats.create('avgPhotEn', _wfr.avgPhotEn)
        ats.create('presCA', _wfr.presCA)
        ats.create('presFT', _wfr.presFT)
        ats.create('unitElFld', _wfr.unitElFld)
        ats.create('unitElFldAng', _wfr.unitElFldAng)

        #Write electric field data
        dataType = np.float32
        if(_form == 1):
            dataType = np.complex64
            
        if(_wfr.arEx is not None): 
            arEx = _wfr.arEx
            if(not isinstance(arEx, (np.ndarray))):
                arEx = np.array(_wfr.arEx)
                
            if(_form == 1):
                arE = np.array([0]*(2*nenxny), dtype=np.float32)
                arEC = arE.view(np.complex64)
                arECR = arEC.reshape(ne, nxny)
                arExC = arEx.view(dataType)
                for i in range(ne):
                    arECR[i] = arExC[i:nenxny:ne]
                arEx = arECR
                    
            data = arEx.view(dataType)
            hf.create_dataset('arEx', data=data, dtype=dataType, shape=data.shape)
            
        if(_wfr.arEy is not None): 
            arEy = _wfr.arEy
            if(not isinstance(arEy, (np.ndarray))):
                arEy = np.array(_wfr.arEy)
 
            if(_form == 1):
                if(arECR is None):
                    arE = np.array([0]*(2*nenxny), dtype=np.float32)
                    arEC = arE.view(np.complex64)
                    arECR = arEC.reshape(ne, nxny)
                    
                arEyC = arEy.view(dataType)
                for i in range(ne):
                    arECR[i] = arEyC[i:nenxny:ne]
                arEy = arECR
                
            data = arEy.view(dataType)
            hf.create_dataset('arEy', data=data, dtype=dataType, shape=data.shape)
            #data = _wfr.arEy.view(np.float32)
            #hf.create_dataset('arEy', data=data, dtype=np.float32, shape=data.shape)

    hf.close()
    return 

#**********************Auxiliary function to write auxiliary/debugging information to an ASCII file:
##def srwl_uti_save_text(_text, _file_path):
##    f = open(_file_path, 'w')
##    f.write(_text + '\n')
##    f.close()
#def srwl_uti_save_text(_text, _file_path, mode='a', newline='\n'): #MR28092016
def srwl_uti_save_text(_text, _file_path, mode='w', newline='\n'): #MR29092016
    with open(_file_path, mode) as f:  
        f.write(_text + newline)

#**********************Auxiliary function to read-in data comumns from ASCII file (2D table):
def srwl_uti_read_data_cols(_file_path, _str_sep, _i_col_start=0, _i_col_end=-1, _n_line_skip=0):
    """
    Auxiliary function to read-in data comumns from ASCII file (2D table)
    :param _file_path: full path (including file name) to the file
    :param _str_sep: column separation symbol(s) (string)
    :param _i_col_start: initial data column to read
    :param _i_col_end: final data column to read
    :param _n_line_skip: number of lines to skip in the beginning of the file
    :return: 2D list containing data columns read
    """
    f = open(_file_path, 'r')
    lines = f.readlines()

    resCols = []

    #nCol = _i_col_end - _i_col_start + 1
    #for iCol in range(nCol):
    #    resCols.append([])

    nRows = len(lines) - _n_line_skip

    for i in range(nRows):
        curLine = lines[_n_line_skip + i]
        #print(curLine)
        
        curLineParts = curLine.split(_str_sep)
        curNumParts = len(curLineParts)
        #print(curLineParts)

        colCount = 0; colCountTrue = 0
        for iCol in range(curNumParts):
            curPart = curLineParts[iCol]
            #print(curPart)
            
            if(len(curPart) > 0):
                if(((_i_col_start <= colCount) or (_i_col_start < 0)) and ((colCount <= _i_col_end) or (_i_col_end < 0))):
                    if len(resCols) < (colCountTrue + 1): resCols.append([])
                    resCols[colCountTrue].append(float(curPart))
                    colCountTrue += 1
                colCount += 1
    f.close()
    return resCols #attn: returns lists, not arrays!

#**********************Auxiliary function to write (save) data comumns to ASCII file (2D table):
def srwl_uti_write_data_cols(_file_path, _cols, _str_sep, _str_head=None, _i_col_start=0, _i_col_end=-1):
    """
    Auxiliary function to write tabulated data (columns, i.e 2D table) to ASCII file
    :param _file_path: full path (including file name) to the file to be (over-)written
    :param _cols: array of data columns to be saves to file
    :param _str_sep: column separation symbol(s) (string)
    :param _str_head: header (string) to write before data columns
    :param _i_col_start: initial data column to write
    :param _i_col_end: final data column to write
    """
    f = open(_file_path, 'w')

    if(_str_head is not None):
        lenStrHead = len(_str_head)
        if(lenStrHead > 0):
            strHead = _str_head
            if(_str_head[lenStrHead - 1] != '\n'):
                strHead = copy(_str_head) + '\n'
            f.write(strHead)
    if(_cols is None):
        f.close(); return
        
    nCols = len(_cols)
    if(nCols <= 0):
        f.close(); return

    nLines = len(_cols[0])
    for i in range(1, nCols):
        newLen = len(_cols[i])
        if(nLines < newLen): nLines = newLen

    strSep = '\t'
    if(_str_sep is not None):
        if(len(_str_sep) > 0): strSep = _str_sep

    strTot = ''
    iColEndP1 = nCols
    if((_i_col_end >= 0) and (_i_col_end < nCols)): iColEndP1 = _i_col_end + 1
    iColEnd = iColEndP1 - 1
    nLinesM1 = nLines - 1
        
    for i in range(nLines):
        curLine = ''
        for j in range(_i_col_start, iColEndP1):
            curElem = ' '
            if(i < len(_cols[j])): curElem = repr(_cols[j][i])
            curLine += curElem
            if(j < iColEnd): curLine += strSep
        if(i < nLinesM1): curLine += '\n'
        strTot += curLine
        
    f.write(strTot)
    f.close()

#**********************Auxiliary function to write auxiliary information about calculation status with "srwl_wfr_emit_prop_multi_e" to an ASCII files:
#def srwl_uti_save_status_mpi(  # #MR20160908
def srwl_uti_save_stat_wfr_emit_prop_multi_e(  #MR20160908
        particle_number=0,
        total_num_of_particles=0,
        #filename='srw_mpi',
        filename='srwl_stat_wfr_emit_prop_multi_e', #MR13012017
        cores=None,
        particles_per_iteration=None
):
    """The function to save .log and .json status files to monitor parallel MPI jobs progress.

    :param particle_number: current particle number.
    :param total_num_of_particles: total number of particles.
    :param filename: the name without extension used to save log/status files.
    :param cores: number of cores used for parallel calculation.
    :param particles_per_iteration: number of particles averaged per iteration (between mpi-receives).-
    :return: None.
    """
    offset = len(str(total_num_of_particles))
    timestamp = '{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
    progress = float(particle_number) / float(total_num_of_particles) * 100.0
    status = 'Running' if progress < 100.0 else 'Finished'

    # Save a log file to monitor duration of calculations:
    mode = 'a'
    if particle_number == 0:
        mode = 'w'
        #text_to_save = '[{}]: Calculation on {} cores with averaging of {} particles/iteration.'.format(
        text_to_save = '[{}]: Calculation on {} core{} with averaging of {} particle{}/iteration.'.format( #OC26042022
            timestamp,
            cores,
            's' if(cores > 1) else '', #OC26042022
            particles_per_iteration,
            's' if(particles_per_iteration > 1) else '', #OC26042022
        )
    else:
        text_to_save = '[{}]: {:8s} {:{offset}d} out of {:{offset}d} ({:6.2f}% complete)'.format(
            timestamp,
            status,
            particle_number,
            total_num_of_particles,
            progress,
            offset=offset,
        )
    status_text_file = '{}.log'.format(filename)
    srwl_uti_save_text(text_to_save, status_text_file, mode=mode)

    # Save JSON file for Sirepo:
    status = {
        'timestamp': timestamp,
        'particle_number': particle_number,
        'total_num_of_particles': total_num_of_particles,
        'progress': progress,
        'status': status,
    }
    status_json_file = '{}.json'.format(filename)
    #with open(status_json_file, 'w') as f:
    #    json.dump(status, f, indent=4, separators=(',', ': '), sort_keys=True)
    #MR05122016: A temp file is created on the same filesystem as the status file to ensure the atomic rename operation:  
    # See the following discussions:  
    # - https://github.com/radiasoft/sirepo/issues/555  
    # - https://github.com/radiasoft/SRW-light/commit/987547f4da3079626bef92aad2f9ca3bade84ada  
    # - http://stackoverflow.com/a/3716361/4143531  
    tmp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', dir=os.path.dirname(status_json_file))
    json.dump(status, tmp_file, indent=4, separators=(',', ': '), sort_keys=True)
    tmp_file.close()
    shutil.move(tmp_file.name, status_json_file)

#**********************Auxiliary function to initialize parameters for the SRW status files and generate the files
def srwl_uti_save_stat_wfr_emit_prop_multi_e_init(num_of_proc, num_part_avg_proc, _tot_num_part, _i_gr=None): #OC31102021 (to be called only by rankMaster!)
#def srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, num_of_proc, num_part_per_proc, num_sent_per_proc, num_part_avg_proc, _tot_num_part=None, _i_gr=None): #OC02032021  
#def srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, num_of_proc, num_part_per_proc, num_sent_per_proc, num_part_avg_proc, _tot_num_part=None): #OC21012021  
#def srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, num_of_proc, num_part_per_proc, num_sent_per_proc, num_part_avg_proc): #MR20012017  
    """Initialize parameters for the SRW status files and generate the files.  
    :param rank: rank of the process.  
    :param num_of_proc: total number of processes.  
    :param num_part_per_proc: number of electrons treated by each worker process.  
    :param num_sent_per_proc: number of sending acts made by each worker process.  
    :param num_part_avg_proc: number of macro-electrons to be used in calculation by each worker.  
    :param _tot_num_part: total number of macro-electrons to be used.  
    :param _i_gr: index of Group of MPI processes.  
    """  
    log_dir = os.path.abspath('__srwl_logs__')  
    #if rank == 0:  
    #OC31102021 (removed the condition above)
    try:  
        os.mkdir(log_dir)  
    except OSError as exc:  
        if exc.errno == errno.EEXIST and os.path.isdir(log_dir):
            pass  
        else:
            raise  
    timestamp = '{:%Y-%m-%d_%H-%M-%S}'.format(datetime.datetime.now()) 
    #log_file = 'srwl_stat_wfr_emit_prop_multi_e_{}'.format(timestamp)
    log_file = 'srwl_stat_wfr_emit_prop_multi_e_' #OC02032021
    if(_i_gr is not None): log_file += repr(_i_gr) #OC02032021
    log_file += '{}' #OC02032021
    log_file = log_file.format(timestamp) #OC02032021

    log_path = os.path.join(log_dir, log_file)

    total_num_of_particles = _tot_num_part #OC31102021 #OC21012021
    #if total_num_of_particles is None: #OC21012021
    #    if num_of_proc <= 1:  
    #        total_num_of_particles = num_part_per_proc  
    #    else:  
    #        total_num_of_particles = num_sent_per_proc * (num_of_proc - 1) * num_part_avg_proc
        
    #if rank == 0:
    #OC31102021 (removed the condition above)
    srwl_uti_save_stat_wfr_emit_prop_multi_e(0, total_num_of_particles, filename=log_path, cores=num_of_proc, particles_per_iteration=num_part_avg_proc)  

    return log_path #OC31102021
    #return log_path, total_num_of_particles

#**********************Auxiliary function to read tabulated 3D Magnetic Field data from ASCII file:
def srwl_uti_read_mag_fld_3d(_fpath, _scom='#'):
    f = open(_fpath, 'r')
    f.readline() #1st line: just pass
    xStart = float(f.readline().split(_scom, 2)[1]) #2nd line: initial X position [m]; it will not actually be used
    xStep = float(f.readline().split(_scom, 2)[1]) #3rd line: step vs X [m]
    xNp = int(f.readline().split(_scom, 2)[1]) #4th line: number of points vs X
    yStart = float(f.readline().split(_scom, 2)[1]) #5th line: initial Y position [m]; it will not actually be used
    yStep = float(f.readline().split(_scom, 2)[1]) #6th line: step vs Y [m]
    yNp = int(f.readline().split(_scom, 2)[1]) #7th line: number of points vs Y
    zStart = float(f.readline().split(_scom, 2)[1]) #8th line: initial Z position [m]; it will not actually be used
    zStep = float(f.readline().split(_scom, 2)[1]) #9th line: step vs Z [m]
    zNp = int(f.readline().split(_scom, 2)[1]) #10th line: number of points vs Z
    totNp = xNp*yNp*zNp
    locArBx = array('d', [0]*totNp)
    locArBy = array('d', [0]*totNp)
    locArBz = array('d', [0]*totNp)
    strSep = '\t'
    for i in range(totNp):
        curLineParts = f.readline().split(strSep)
        #DEBUG
        #print(i, curLineParts)
        #END DEBUG
        locArBx[i] = float(curLineParts[0])
        locArBy[i] = float(curLineParts[1])
        locArBz[i] = float(curLineParts[2])
    f.close()
    xRange = xStep
    if xNp > 1: xRange = (xNp - 1)*xStep
    yRange = yStep
    if yNp > 1: yRange = (yNp - 1)*yStep
    zRange = zStep
    if zNp > 1: zRange = (zNp - 1)*zStep
    
    xc = xStart + 0.5*xStep*(xNp - 1)
    yc = yStart + 0.5*yStep*(yNp - 1)
    zc = zStart + 0.5*zStep*(zNp - 1)
    return SRWLMagFldC(SRWLMagFld3D(locArBx, locArBy, locArBz, xNp, yNp, zNp, xRange, yRange, zRange, 1), xc, yc, zc)

#**********************Auxiliary function to allocate array
#(to walk-around the problem that simple allocation "array(type, [0]*n)" at large n is usually very time-consuming)
def srwl_uti_array_alloc(_type, _n, _list_base=[0]): #OC14042019
#def srwl_uti_array_alloc(_type, _n):

    # import numpy as np #OCTEST OC02112022
    # resAr = np.zeros(_n, dtype=float) #OCTEST OC02112022
    # resAr = array(_type, _list_base*_n) #OCTEST OC02112022

    #nPartMax = 1000000 #OCTEST OC02112022
    nPartMax = 10000000 #to tune
    #nPartMax = 3000000 #to tune
    
    #print('srwl_uti_array_alloc: array requested:', _n)
    lenBase = len(_list_base) #OC14042019
    nTrue = _n*lenBase #OC14042019
    
    if(nTrue <= nPartMax): return array(_type, _list_base*_n)
    #if(_n <= nPartMax): return array(_type, [0]*_n)
        #resAr = array(_type, [0]*_n)
        #print('Array requested:', _n, 'Allocated:', len(resAr))
        #return resAr

    nPartMax_d_lenBase = int(nPartMax/lenBase) #OC14042019

    nEqualParts = int(_n/nPartMax_d_lenBase) #OC14042019
    #nEqualParts = int(_n/nPartMax)

    nResid = int(_n - nEqualParts*nPartMax_d_lenBase) #OC14042019
    #nResid = int(_n - nEqualParts*nPartMax)

    resAr = array(_type, _list_base*nPartMax_d_lenBase) #OC14042019
    #resAr = array(_type, [0]*nPartMax)

    if(nEqualParts > 1):
        auxAr = deepcopy(resAr)
        for i in range(nEqualParts - 1): 
            resAr.extend(auxAr)
    if(nResid > 0):
        auxAr = array(_type, _list_base*nResid) #OC14042019
        #auxAr = array(_type, [0]*nResid)
        resAr.extend(auxAr)

    #print('Array requested:', _n, 'Allocated:', len(resAr))
    
    return resAr

#**********************Auxiliary function to generate Halton sequence (to replace pseudo-random numbers)
#Contribution from R. Lindberg, X. Shi (APS)
def srwl_uti_math_seq_halton(i, base=2):
#def Halton(i, base=2):
    h = 0
    fac = 1.0/base
    while i != 0:
        digit = i % base
        #h = h + digit*fac
        h += digit*fac
        i = (i - digit)/base
        #fac = fac/base
        fac /= base
    return h

#****************************************************************************
#****************************************************************************
#Wavefront manipulation functions
#****************************************************************************
#****************************************************************************

#**********************Auxiliary function to propagate wavefront over free space in a number of steps and extract resulting intensity cuts
def srwl_wfr_prop_drifts(_wfr, _dz, _nz, _pp, _do3d=False, _nx=-1, _ny=-1, _rx=0, _ry=0, _xc=0, _yc=0, _pol=6, _type=0, _ord_interp=1):
    """
    Propagates wavefront over free space in a number of steps and generates intensity distributions vs (z,x) and (z,y) and possibly (z,x,y)
    :param _wfr: input/output wavefront
    :param _dz: longitudinal step size for the drifts
    :param _nz: number of drift steps to be made
    :param _pp: list of propagation parameters to be used for each drift step
    :param _do3d: generate intensity vs (z,x,y) in addition to the cuts (z,x) and (z,y)
    :param _nx: number of points vs horizontal position (is taken into account if >0)
    :param _ny: number of points vs vertical position (is taken into account if >0)
    :param _rx: number of points vs horizontal position (is taken into account if >0)
    :param _ry: number of points vs vertical position(is taken into account if >0)
    :param _xc: horizontal position for vertical cut of intensity
    :param _yc: vertical position for horizontal cut of intensity
    :param _pol: switch specifying polarization component to be extracted:
            =0 -Linear Horizontal; 
            =1 -Linear Vertical; 
            =2 -Linear 45 degrees; 
            =3 -Linear 135 degrees;
            =4 -Circular Right; 
            =5 -Circular Left; 
            =6 -Total
    :param _type: switch specifying "type" of a characteristic to be extracted:
            =0 -"Single-Electron" Intensity; 
            =1 -"Multi-Electron" Intensity; 
            =2 -"Single-Electron" Flux; 
            =3 -"Multi-Electron" Flux; 
            =4 -"Single-Electron" Radiation Phase; 
            =5 -Re(E): Real part of Single-Electron Electric Field;
            =6 -Im(E): Imaginary part of Single-Electron Electric Field;
            =7 -"Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence)
    :param _ord_interp: interpolation order for final intensity calc.
    :return resulting intensity distributions
    """

    if(_nx <= 0): _nx = _wfr.mesh.nx
    if(_ny <= 0): _ny = _wfr.mesh.ny

    nzp1 = _nz + 1

    resIntVsZX = array('f', [0]*(nzp1*_nx)) #array to store the final intensity
    resIntVsZY = array('f', [0]*(nzp1*_ny)) #array to store the final intensity
    resIntVsZXY = None
    if(_do3d): resIntVsZXY = array('f', [0]*(nzp1*_nx*_ny)) #array to store the final intensity

    #OC19102018
    resFWHMxVsZ = array('f', [0]*nzp1)
    resFWHMyVsZ = array('f', [0]*nzp1)

    ec = _wfr.mesh.eStart
    if(_wfr.mesh.ne > 1): ec = 0.5*(_wfr.mesh.eStart + _wfr.mesh.eFin)

    cntDrift = SRWLOptC([SRWLOptD(_dz)], [_pp])

    xStart = _wfr.mesh.xStart
    xFin = _wfr.mesh.xFin
    if(_rx > 0):
        halfRx = 0.5*_rx
        xStart = _xc - halfRx
        xFin = _xc + halfRx
    xStep = (xFin - xStart)/(_nx - 1) if _nx > 1 else 0

    yStart = _wfr.mesh.yStart
    yFin = _wfr.mesh.yFin
    if(_ry > 0):
        halfRy = 0.5*_ry
        yStart = _yc - halfRy
        yFin = _yc + halfRy
    yStep = (yFin - yStart)/(_ny - 1) if _ny > 1 else 0

    #OC19102018
    meshFWHMx = copy(_wfr.mesh)
    meshFWHMx.ne = 1; meshFWHMx.ny = 1
    meshFWHMx.eStart = ec; meshFWHMx.eFin = ec
    meshFWHMx.yStart = _yc; meshFWHMx.yFin = _yc
    meshFWHMy = copy(_wfr.mesh)
    meshFWHMy.ne = 1; meshFWHMy.nx = 1
    meshFWHMy.eStart = ec; meshFWHMy.eFin = ec
    meshFWHMy.xStart = _xc; meshFWHMy.xFin = _xc

    for iz in range(0, nzp1):

        if(iz > 0):
            print('Propagation (step # ' + repr(iz) + ') ... ', end='')
            t0 = time.time();
            srwl.PropagElecField(_wfr, cntDrift)
            print('completed (lasted', round(time.time() - t0, 6), 's)')

        curMesh = _wfr.mesh
        xStepCurMesh = (curMesh.xFin - curMesh.xStart)/(curMesh.nx - 1)
        yStepCurMesh = (curMesh.yFin - curMesh.yStart)/(curMesh.ny - 1)
        curIntVsX = array('f', [0]*curMesh.nx)
        curIntVsY = array('f', [0]*curMesh.ny)
        
        srwl.CalcIntFromElecField(curIntVsX, _wfr, _pol, _type, 1, ec, _xc, _yc) #int. vs X

        #OC19102018
        meshFWHMx.nx = curMesh.nx; meshFWHMx.xStart = curMesh.xStart; meshFWHMx.xFin = curMesh.xFin
        intInfX = srwl.UtiIntInf(curIntVsX, meshFWHMx)
        #print(intInfX)
        resFWHMxVsZ[iz] = intInfX[4]

        x = xStart
        for ix in range(_nx): #interpolation
            resIntVsZX[iz + ix*nzp1] = uti_math.interp_1d(x, curMesh.xStart, xStepCurMesh, curMesh.nx, curIntVsX, _ord_interp)
            x += xStep
        
        srwl.CalcIntFromElecField(curIntVsY, _wfr, _pol, _type, 2, ec, _xc, _yc) #int. vs Y

        #OC19102018
        meshFWHMy.ny = curMesh.ny; meshFWHMy.yStart = curMesh.yStart; meshFWHMy.yFin = curMesh.yFin
        intInfY = srwl.UtiIntInf(curIntVsY, meshFWHMy)
        #print(intInfY)
        resFWHMyVsZ[iz] = intInfY[4]

        y = yStart
        for iy in range(_ny): #interpolation
            resIntVsZY[iz + iy*nzp1] = uti_math.interp_1d(y, curMesh.yStart, yStepCurMesh, curMesh.ny, curIntVsY, _ord_interp)
            y += yStep

        if(_do3d):
            curIntVsXY = array('f', [0]*_wfr.mesh.nx*_wfr.mesh.ny)
            srwl.CalcIntFromElecField(curIntVsXY, _wfr, _pol, _type, 3, ec, _xc, _yc) #int. vs XY

            nxny = _nx*_ny
            y = yStart
            for iy in range(_ny):
                x = xStart
                for ix in range(_nx):
                    resIntVsZXY[ix + iy*_nx + iz*nxny] = uti_math.interp_2d(x, y, curMesh.xStart, xStepCurMesh, curMesh.nx, curMesh.yStart, yStepCurMesh, curMesh.ny, curIntVsXY, _ord_interp)
                    x += xStep
                y += yStep

    resMesh = SRWLRadMesh(ec, ec, 1, xStart, xFin, _nx, yStart, yFin, _ny, 0.)
    resMesh.zFin = _dz*_nz #adding long. mesh params (note: these may not propagate further!)
    resMesh.nz = nzp1 #adding long. mesh params (may not propagate further!)

    #return resIntVsZX, resIntVsZY, resIntVsZXY, resMesh
    return resIntVsZX, resIntVsZY, resIntVsZXY, resMesh, resFWHMxVsZ, resFWHMyVsZ #OC19102018

#**********************Auxiliary function to setup Coherent Wavefront from Intensity (assuming spherical or astigmatic wave)
def srwl_wfr_from_intens(_ar_int, _mesh, _part_beam, _Rx, _Ry, _xc=0, _yc=0):
    """
    Setup Coherent Wavefront from Intensity (assuming spherical or asigmatic wave): note this is an error-prone procedure
    :param _ar_int: input intensity array 
    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
    :param _part_beam: Finite-Emittance beam (SRWLPartBeam type)
    :param _Rx: horizontal wavefront radius [m]
    :param _Ry: vertical wavefront radius [m]
    :param _xc: horizontal wavefront center position [m]
    :param _yc: vertical wavefront center position [m]
    """

    lenInt = len(_ar_int)
    nTot = _mesh.ne*_mesh.nx*_mesh.ny
    if(lenInt != nTot):
        raise Exception("Mesh parameters are not consistent with the length of intensity array") 
    
    aux_const = 3.14159265358979E+06/1.23984186
    constRx = aux_const/_Rx
    constRy = aux_const/_Ry

    wfr = SRWLWfr()
    wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions (may be modified by the library!)

    eStep = 0 if(_mesh.ne <= 1) else (_mesh.eFin - _mesh.eStart)/(_mesh.ne - 1)
    xStep = 0 if(_mesh.nx <= 1) else (_mesh.xFin - _mesh.xStart)/(_mesh.nx - 1)
    yStep = 0 if(_mesh.ny <= 1) else (_mesh.yFin - _mesh.yStart)/(_mesh.ny - 1)

    halfPerX = _mesh.ne
    halfPerY = halfPerX*_mesh.nx
    ePh = _mesh.eStart

    for ie in range(_mesh.ne):
        constRxE = constRx*ePh
        constRyE = constRy*ePh

        y = _mesh.yStart - _yc
        for iy in range(_mesh.ny):
            dPhaseY = constRyE*y*y

            halfPerYiy_p_ie = halfPerY*iy + ie

            x = _mesh.xStart - _xc
            for ix in range(_mesh.nx):
                phase = dPhaseY + constRxE*x*x
                cosPhase = cos(phase)
                sinPhase = sin(phase)
                
                ofstI = halfPerYiy_p_ie + halfPerX*ix
                curI = abs(_ar_int[ofstI])
                magn = sqrt(curI)
                
                ofstE = ofstI*2
                wfr.arEx[ofstE] = magn*cosPhase
                wfr.arEy[ofstE] = 0
                ofstE += 1
                wfr.arEx[ofstE] = magn*sinPhase
                wfr.arEy[ofstE] = 0

                x += xStep
            y += yStep
        ePh += eStep

    #More wavefront parameters
    wfr.partBeam = _part_beam
    wfr.mesh = deepcopy(_mesh)
    wfr.Rx = _Rx #instant wavefront radii
    wfr.Ry = _Ry
    wfr.dRx = 0.01*_Rx #error of wavefront radii
    wfr.dRy = 0.01*_Ry
    wfr.xc = _xc #instant transverse coordinates of wavefront instant "source center"
    wfr.yc = _yc
    wfr.avgPhotEn = _mesh.eStart if(_mesh.ne == 1) else 0.5*(_mesh.eStart + _mesh.eFin) #average photon energy for time-domain simulations
    wfr.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
    wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
    wfr.unitElFld = 1 #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)
    #wfr.arElecPropMatr = array('d', [0] * 20) #effective 1st order "propagation matrix" for electron beam parameters
    #wfr.arMomX = array('d', [0] * 11 * _ne) #statistical moments (of Wigner distribution); to check the exact number of moments required
    #wfr.arMomY = array('d', [0] * 11 * _ne)
    #wfr.arWfrAuxData = array('d', [0] * 30) #array of auxiliary wavefront data

    return wfr

#**********************Auxiliary function to generate standard filenames for radiation characteristics
def srwl_wfr_fn(_fn_core, _type, _form='ascii'): #OC15022021
#def srwl_wfr_fn(_fn_core, _type):
    """
    Generate standard filenames for radiation characteristics from a given core name
    :param _fn_core: core filename
    :param _type: type of radiation characteristic:
        0- Electric Field
        1- Spectral Intensity, i.e. Spectral Flux per Unit Surface Area
        2- Spectral Angular Intensity, i.e. Spectral Flux per Unit Solid Angle
        3- Mutual Intensity (in real space, i.e. dependent on coordinates)
        31- Mutual Intensity Cut vs X
        32- Mutual Intensity Cut vs Y
        4- Degree of Coherence
        41- Degree of Coherence Cut vs X
        42- Degree of Coherence Cut vs Y
        5- Wigner Distribution / Brightness
        51- Wigner Distribution / Brightness Cut vs X
        52- Wigner Distribution / Brightness Cut vs Y
        7- Coherent Modes
    :param _form: format of data file ('ascii' and 'hdf5' or 'h5' are supported)
    """

    if(_type < 0): return _fn_core 

    #stdExt = '.dat'
    #OC15022021
    stdExt = '.dat'
    if((_form == 'hdf5') or (_form == 'h5')): stdExt = '.h5'
    fnCore = copy(_fn_core)

    len_fnCore = len(fnCore)
    #ext_fnCore = fnCore[(len_fnCore - 4):len_fnCore]
    #if((ext_fnCore == '.txt') or (ext_fnCore == '.dat')):
    #    fnCore = fnCore[0:(len_fnCore - 4)]
    #    stdExt = ext_fnCore
    #OC14022021
    indDot = fnCore.rfind('.')

    #DEBUG
    #print(fnCore, indDot, len_fnCore)
    #END DEBUG
    
    if(indDot >= 0):
        fnCore = fnCore[0:indDot]
        stdExt = _fn_core[indDot:len_fnCore]

        #DEBUG
        #print(fnCore, stdExt)
        #END DEBUG

    suf = ''
    if(_type == 0): suf = '_ef'
    elif(_type == 1): suf = '_ic'
    elif(_type == 2): suf = '_ia'
    elif(_type == 3): suf = '_mi'
    elif(_type == 31): suf = '_mix'
    elif(_type == 32): suf = '_miy'
    elif(_type == 4): suf = '_dc'
    elif(_type == 41): suf = '_dcx'
    elif(_type == 42): suf = '_dcy'
    elif(_type == 5): suf = '_wd'
    elif(_type == 51): suf = '_wdx'
    elif(_type == 52): suf = '_wdy'
    elif(_type == 7): suf = '_cm'

    return fnCore + suf + stdExt

#**********************Auxiliary function to check/correct/change file extension, depending on requested file format and type of calculation
def srwl_wfr_fn_ext(_fp, _form, _char=0): #OC19072021
    """
    Check / change file extension, depending on requested file format and type of calculation
    :param _fp: filename (/path)
    :param _form: required file format, can be:
        'ascii' or 'ASCII' or 'asc' or 'ASC' for ASCII format
        'hdf5' or 'HDF5' or 'h5' or 'H5' for HDF5 format
    :param _char: radiation characteristic to be calculated (one of values supported by the srwl_wfr_emit_prop_multi_e function)
    :return: updated / checked filename
    """

    if((_fp is None) or (_form is None)): return _fp

    reqForm = copy(_form)
    reqForm = reqForm.upper()
    if(reqForm == 'ASCII'): reqForm = 'ASC'
    elif(reqForm == 'HDF5'): reqForm = 'H5'
    
    if(_char in [6, 61, 7]): reqForm = 'H5' #for the above cases, only HDF5 format is currently supported

    reqExt = None
    if(reqForm == 'H5'): reqExt = '.h5'
    elif(reqForm == 'ASC'): reqExt = '.dat'
    else: reqExt = '' #??

    lenFP = len(_fp)
    indDot = _fp.rfind('.')
    fnCore = copy(_fp)
    origExt = ''
    if(indDot >= 0):
        fnCore = _fp[0:indDot]
        origExt = _fp[indDot:lenFP]

    if(origExt == reqExt): return _fp
    else: return fnCore + reqExt

#**********************Auxiliary function to post-process (average) CSD data from files generated e.g. by different groups of MPI processes
def srwl_wfr_csd_avg(_fp_core, _fp_ext, _fi_st, _fi_en, _csd0=None, _awfr=None, _form='hdf5', _do_fin_save=True, _do_del_aux_files=False, _do_sym=True, _fin_fp_core=None, _do_sum=False): #OC27102021
#def srwl_wfr_csd_avg(_fp_core, _fp_ext, _fi_st, _fi_en, _csd0=None, _awfr=None, _form='hdf5', _do_fin_save=True, _do_del_aux_files=False, _do_sym=True, _fin_fp_core=None): #OC11102021
#def srwl_wfr_csd_avg(_fp_core, _fp_ext, _fi_st, _fi_en, _csd0=None, _awfr=None, _form='hdf5', _do_fin_save=True, _do_del_aux_files=False, _do_sym=True): #OC21062021
#def srwl_wfr_csd_avg(_fp_core, _fp_ext, _fi_st, _fi_en, _csd0=None, _awfr=None, _form='hdf5', _do_fin_save=True, _do_del_aux_files=False): #OC20062021
#def srwl_wfr_csd_avg(_fp_core, _fp_ext, _fi_st, _fi_en, _csd0=None, _awfr=None, _form='hdf5', _do_fin_save=True): #OC19062021

    arCSD = None
    meshCSD = None
    fi_st = _fi_st
    if(_csd0 is not None):
        arCSD = _csd0.arS
        meshCSD = _csd0.mesh
        fi_st = _fi_st + 1

    awfr = _awfr #OC10102021
    for j in range(fi_st, _fi_en+1):

        curFP = _fp_core + '_' + repr(j) + _fp_ext
        #curFP = _fp_core + repr(j) + _fp_ext
        curFP = srwl_wfr_fn(curFP, 3) #adds suffix "_mi" before extension

        #DEBUG
        t0 = time.time()
        #END DEBUG

        resPartMI = srwl_uti_read_intens(_file_path=curFP, _form=_form)
        arPartMI = resPartMI[0]
        meshPartMI = resPartMI[1] #The mesh is supposed to be the same for all partial CSDs
        extPar = resPartMI[2]
        if(awfr is None): #OC10102021
            if(j == fi_st): awfr = resPartMI[3]

        #DEBUG
        t1 = round(time.time() - t0)
        print('Reading-in of partial CSD file #:', j, 'accomplished in:', t1, 's')
        sys.stdout.flush()
        #END DEBUG

        if(j == _fi_st):
            arCSD = arPartMI
            meshCSD = meshPartMI
        
        else: #OC04102021

            #DEBUG
            #print('Adding CSD data from Group:', j, 'is about to start')
            t0 = time.time()
            #sys.stdout.flush()
            #END DEBUG

            itNum = j - fi_st + 1
            if _do_sum: itNum = -1 #OC27102021
            
            srwl.UtiIntProc(arCSD, meshCSD, arPartMI, meshPartMI, [1, itNum]) #Averaging or Summing-up CSDs
            #srwl.UtiIntProc(arCSD, meshCSD, arPartMI, meshPartMI, [1, (j - fi_st + 1)]) #Averaging CSDs
            #srwl.UtiIntProc(arCSD, meshCSD, arPartMI, meshPartMI, [1, j]) #Averaging CSDs

            #DEBUG
            print('Partial CSD data from MPI Group #', j, 'was added in:', round(time.time() - t0), 's')
            sys.stdout.flush()
            #END DEBUG

            if(_do_del_aux_files): os.remove(curFP)

        if(j > _fi_st): del resPartMI

    if(_do_sym): 

        #DEBUG
        t0 = time.time()
        #END DEBUG

        srwl.UtiIntProc(arCSD, meshCSD, None, None, [4]) #Fill-in "symmetrical" part of the Hermitian MI distribution

        #DEBUG
        print('Symmetrical parts of the CSD data (Hermitian matrix) were filled-out in:', round(time.time() - t0), 's')
        sys.stdout.flush()
        #END DEBUG

    if(_do_fin_save): #Final saving is done here

        fpResMI = srwl_wfr_fn(_fp_core + _fp_ext, 3) #adds suffix "_mi" before extension
        #fpResMI = _fp_core + _fp_ext
        if(_fin_fp_core is not None): #OC11102021
            if(len(_fin_fp_core) > 0):
                fpResMI = srwl_wfr_fn(_fin_fp_core + _fp_ext, 3) #adds suffix "_mi" before extension

        numComp = 1
        arLabels = ['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity']
        arUnits = ['eV', 'm', 'm', 'ph/s/.1%bw/mm^2']
        mutual = 1
        cmplx = 1
        if(extPar is not None):
            numComp = extPar['n_stokes']
            arLabels = extPar['arLabels']
            arUnits = extPar['arUnits']
            mutual = extPar['mutual']
            cmplx = extPar['cmplx']

        #DEBUG
        t0 = time.time()
        #END DEBUG

        srwl_uti_save_intens(arCSD, meshCSD, fpResMI, numComp, _arLabels = arLabels, _arUnits = arUnits, _mutual = mutual, _cmplx = cmplx, _form = _form, _wfr = awfr) #OC10102021
        #srwl_uti_save_intens(arCSD, meshCSD, fpResMI, numComp, _arLabels = arLabels, _arUnits = arUnits, _mutual = mutual, _cmplx = cmplx, _form = _form, _wfr = _awfr)

        #DEBUG
        print('Total CSD data saved in:', round(time.time() - t0), 's')
        sys.stdout.flush()
        #END DEBUG

    if(_csd0 is not None): return _csd0, awfr #OC10102021
    #if(_csd0 is not None): return _csd0
    else: 
        csd = SRWLStokes(_arS=arCSD, _typeStokes='f', _mutual=1, _n_comp=1)
        csd.mesh = meshCSD
        return csd, awfr #OC10102021
        #return csd

#**********************Auxiliary function to perform CMD on 4D CSD data
def srwl_wfr_cmd(_csd, _n_modes, _awfr=None, _alg=None): #OC21062021
    """
    Perform Coherent Mode Decomposition (CMD) of (4D) Cross-Spectral Density / Mutual Intensity data.
    :param _csd: Cross-Spectral Density / Mutual Intensity described by SRWLStokes object with complex 4D data array
    :param _n_modes: maximum number of coherent modes to be produced by the decomposition
    :param _awfr: optional average Wavefront that can be used if quadratic phase terms have to be treated during the decomposition
    :param _alg: algorithm to be used for solving the CMD Eigenvalue problem; among those currently supported are: 'PM' - Primme (requires a dedicated 'primme' module to be installed); 'SP' - SciPy; 'SPS' - SciPy Sparse (the latter two available in SciPy package)
    :return: list of Coherent Modes (objects of SRWLWfr type) produced by the decomposition
    """

    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')

    if(not isinstance(_csd, SRWLStokes)):
        raise Exception('Incorrect CSD object submitted (SRWLStokes type expected)')

    ndd = _csd.mesh.nx*_csd.mesh.ny
    nTot = ndd*ndd*2
    lenS = len(_csd.arS)

    data = None
    if(isinstance(_csd.arS, array) or (lenS > nTot)): #Conversion to NumPy array is required (to check)?

        data = np.array([0]*nTot, 'f')
        data[:] = _csd.arS
        #data[0:nTot] = _csd.arS #?

    elif(isinstance(_csd.arS, np.ndarray)): data = _csd.arS #OC29062021
    #elif(isinstance(_csd.arS, np.array)): data = _csd.arS
    else:
        raise Exception('CSD data should be in Python or NumPy array')

    if(len(_csd.arS) > nTot):
        print('WARNING: CMD will be performed on only one Stokes component')

    treatComQuadPhTerms = False
    if(_awfr is not None):
        if((abs(_awfr.Rx) > 5*_awfr.dRx) and (abs(_awfr.Ry) > 5*_awfr.dRy)): treatComQuadPhTerms = True
        #if((abs(_awfr.Rx) > _awfr.dRx) and (abs(_awfr.Ry) > _awfr.dRy)): treatComQuadPhTerms = True

    if(treatComQuadPhTerms): 
        #DEBUG
        t0 = time.time()
        #END DEBUG

        srwl.UtiIntProc(data, _csd.mesh, None, None, [5, -1, _awfr.Rx, _awfr.Ry, _awfr.xc, _awfr.yc]) #Subtract common quadratic phase terms from CSD

        #DEBUG
        t1 = round(time.time() - t0, 3)
        print('Common Quadratic Terms were subtracted from CSD in:', t1, 's')
        sys.stdout.flush()
        #END DEBUG
       
    #DEBUG
    t0 = time.time()
    #END DEBUG

    data = data.view(np.complex64).reshape(ndd, ndd)

    eigSys = UtiMathEigen(dim=ndd) if(_alg is None) else UtiMathEigen(dim=ndd, alg=_alg)  #OC01072021 (because UtiMathEigen ctor has default alg='PM')

    cohModes, eigVals = eigSys.eigen_left(data, n_modes=_n_modes) #OC29062021
    #cohModes, eigVals = eigSys.eigen_left(data, n_modes=_n_modes, alg=_alg)

    #DEBUG
    t1 = round(time.time() - t0, 3)
    print('Coherent Mode Decomposition done in:', t1, 's')
    sys.stdout.flush()
    #END DEBUG
    
    if(treatComQuadPhTerms):
        #DEBUG
        t0 = time.time()
        #END DEBUG

        srwl.UtiIntProc(cohModes, _csd.mesh, None, None, [6, _n_modes, 1, _awfr.Rx, _awfr.Ry, _awfr.xc, _awfr.yc]) #Add back common quadratic phase terms to CMs

        #DEBUG
        t1 = round(time.time() - t0, 3)
        print('Common Quadratic Terms were added to CMs in:', t1, 's')
        sys.stdout.flush()
        #END DEBUG

    return cohModes, eigVals

#**********************Main Partially-Coherent Emission and Propagaiton simulation function
def srwl_wfr_emit_prop_multi_e(_e_beam, _mag, _mesh, _sr_meth, _sr_rel_prec, _n_part_tot, _n_part_avg_proc=1, _n_save_per=100,
                               _file_path=None, _sr_samp_fact=-1, _opt_bl=None, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0,
                               #_rand_meth=1, _tryToUseMPI=True):
                               #_rand_meth=1, _tryToUseMPI=True, _w_wr=0.): #OC26032016 (added _w_wr)
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0.): #OC07092016 (renamed _wr)
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _det=None): #OC06122016
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None): #OC05012017
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0): #OC05042017
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False): #OC14082018
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False): #OC24042020
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False, _file_form='ascii'): #OC05022021
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False, _file_form='ascii', _com_mpi=None, _n_mpi=1): #OC01032021
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False, _file_form='ascii', _n_mpi=1): #OC02032021
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False, _file_form='ascii', _n_mpi=1, _del_aux_files=False): #OC02032021
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False, _file_form='ascii', _n_mpi=1, _n_cm=1000, _del_aux_files=False): #OC27062021
                               _rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0, _file_bkp=False, _rand_opt=False, _file_form='ascii', _n_mpi=1, _n_cm=1000, _ms=0, _del_aux_files=False): #OC22112022
    """
    Calculate Stokes Parameters of Emitted (and Propagated, if beamline is defined) Partially-Coherent SR.
    :param _e_beam: Finite-Emittance e-beam (SRWLPartBeam type)
    :param _mag: Magnetic Field container (magFldCnt type) or Coherent Gaussian beam (SRWLGsnBm type) or Point Source (SRWLPtSrc type) or list of Coherent Modes (of SRWLWfr type each)
    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
    :param _sr_meth: SR Electric Field calculation method to be used: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    :param _sr_rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, the smaller the more accurate)
    :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
    :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
    :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
    :param _file_path: path to file for saving intermediate average Stokes data by master process
    :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
    :param _opt_bl: optical beamline (container) to propagate the radiation through (SRWLOptC type)
    :param _pres_ang: switch specifying presentation of the resulting Stokes parameters: coordinate (0) or angular (1) or both (2, to implement !!!)
    :param _char: radiation characteristic to calculate:
        0- Total Intensity, i.e. Flux per Unit Surface Area (s0);
        1- Four Stokes components of Flux per Unit Surface Area;
        2- Mutual Intensity Cut (2D) vs X;
        3- Mutual Intensity Cut (2D) vs Y;
        4- Mutual Intensity and Degree of Coherence Cuts (2D) vs X & Y;
        5- Degree of Coherence Cuts (2D) vs X & Y;
        6- Mutual Intensity / Cross-Spectral Density (4D) vs X & Y;
        61- Coherent Modes (array of 2D distributions of Electric Field vs X & Y) and Mutual Intensity / Cross-Spectral Density (4D) vs X & Y it is based on;
        7- Coherent Modes (array of 2D distributions of Electric Field vs X & Y);
        10- Flux
        11- Four Stokes components of Flux (?)
        20- Electric Field (sum of fields from all macro-electrons, assuming CSR)
        40- Total Intensity, i.e. Flux per Unit Surface Area (s0), Mutual Intensity and Degree of Coherence Cuts vs X & Y;
        41- Total Intensity, i.e. Flux per Unit Surface Area (s0) and Degree of Coherence Cuts vs X & Y;
    :param _x0: horizontal center position for mutual intensity calculation
    :param _y0: vertical center position for mutual intensity calculation
    :param _e_ph_integ: integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from _mesh
    :param _rand_meth: method for generation of pseudo-random numbers for e-beam phase-space integration:
        1- standard pseudo-random number generator
        2- Halton sequences
        3- LPtau sequences (to be implemented)
    :param _tryToUseMPI: switch specifying whether MPI should be attempted to be used
    :param _wr: initial wavefront radius [m] to assume at wavefront propagation (is taken into account if != 0)
    :param _wre: initial wavefront radius error [m] to assume at wavefront propagation (is taken into account if != 0)
    :param _det: detector object for post-processing of final intensity (instance of SRWLDet)
    :param _me_approx: approximation to be used at multi-electron integration: 0- none (i.e. do standard M-C integration over 5D phase space volume of e-beam), 1- integrate numerically only over e-beam energy spread and use convolution to treat transverse emittance
    :param _file_bkp: create or not backup files with resulting multi-electron radiation characteristics
    :param _rand_opt: randomize parameters of optical elements at each fully-coherent wavefront propagation (e.g. to simulate impact of vibrations) or not
    :param _file_form: format of output files ('ascii' / 'asc' and 'hdf5' supported)
    :param _n_mpi: number of independent "groups" of MPI processes (for 4D CSD calculation)
    :param _n_cm: number of coherent modes to calculate (is taken into account if _char==61 or _char==7)
    :param _ms: index of the first wavefront / coherent mode to start computation
    :param _del_aux_files: delete (or not) auxiliary files (applies to different types of calculations)
   """

    doMutual = 0 #OC30052017
    #if((_char >= 2) and (_char <= 4)): doMutual = 1
    #if((_char >= 2) and (_char <= 5)): doMutual = 1 #OC15072019
    #if((_char >= 2) and (_char <= 6)): doMutual = 1 #OC02022021
    if(((_char >= 2) and (_char <= 6)) or (_char == 61) or (_char == 7)): doMutual = 1 #OC27062021

    if((_det is not None) and ((doMutual > 0) and ((_char != 6) and (_char != 61) and (_char != 7)))): raise Exception("Detector processing is not supported for mutual intensity") #OC20062021
    #if((_det is not None) and ((doMutual > 0) and (_char != 6))): raise Exception("Detector processing is not supported for mutual intensity") #OC02022021
    #if((_det is not None) and (doMutual > 0)): raise Exception("Detector processing is not supported for mutual intensity") #OC30052017

    #DEBUG
    #print('_mesh.xStart=', _mesh.xStart, '_mesh.xFin=', _mesh.xFin, '_mesh.yStart=', _mesh.yStart, '_mesh.yFin=', _mesh.yFin) #DEBUG
    #print('_sr_samp_fact:', _sr_samp_fact) #DEBUG
    #for i in range(len(_opt_bl.arProp)): #DEBUG
    #    print(i, _opt_bl.arProp[i]) #DEBUG

    #DEBUG
    #self.arOpt = []
    #self.arProp = []
    #print('_file_bkp=',_file_bkp)
    #print('srwl_wfr_emit_prop_multi_e: _e_ph_integ=', _e_ph_integ)

    nProc = 1
    rank = 0 #OC30122021
    #rank = 1
    MPI = None
    comMPI = None

    if(_tryToUseMPI):
        try:
            from mpi4py import MPI #OC091014
            #DEBUG
            #print('mpi4py loaded')
            #END DEBUG

            #if(_com_mpi is None): #OC01032021

            #DEBUG
            #resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], -1) #MPI module load
            #resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], 0) #MPI module load
            #print('__import__ passed')
            #MPI = resImpMPI4Py.MPI
         
            comMPI = MPI.COMM_WORLD

            #else:
            #    comMPI = _com_mpi

            rank = comMPI.Get_rank()
            nProc = comMPI.Get_size()

        except:
            print('Calculation will be sequential (non-parallel), because "mpi4py" module can not be loaded or used') #OC01032021
            #print('Calculation will be sequential (non-parallel), because "mpi4py" module can not be loaded')
    
    nProcTot = nProc #OC30102021
    if(nProcTot > 1): #OC30122021
        if(nProcTot < _n_mpi*2): raise Exception("Number of independent \"groups\" of MPI processes is too large for the requested total number of MPI processes.")
    #if(nProcTot < _n_mpi*2): raise Exception("Number of independent \"groups\" of MPI processes is too large for the requested total number of MPI processes.")

    #DEBUG (use this to mimic MPI-execution)
    #nProc = 4 #8 #41
    #rank = 0 #3 #5 #0 #3 #2 #4 #0 #5
    #END DEBUG
    #print('DEBUG: rank, nProc:', rank, nProc)

    #OC02032021: nProc, rankMaster, itStartEnd, _file_path are modified below; rank is not modified
    rankMaster = 0 #OC02032021: rank of Master
    itStartEnd = None #OC02032021: #start and end indexes for aligned conjugated coordinate
    iGr = 0 #OC26042022 #index of Group (out of _n_mpi), used only at _char = 6
    #iGr = None #index of Group (out of _n_mpi), used only at _char = 6
    fpCore = None
    fpExt = None
    fp_cm = None #OC28062021

    #OC19072021: check/correct/change file extension, depending on requested file format and type of calculation
    _file_path = srwl_wfr_fn_ext(_file_path, _file_form, _char)

    nPartTotOrig = _n_part_tot #OC30102021
    arRankMasters = [0] #OC27042022 (?)
    #arRankMasters = None #OC30102021

    if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021
    #if(_char == 6): 

        fp_cm = srwl_wfr_fn(_file_path, _type=7, _form='hdf5') #OC02102021

        _n_part_avg_proc = 1 #OC18022021

        if((_n_mpi > 1) and (nProc <= 1)): _n_mpi = 1 #OC18062021

        if((_n_mpi > 1) and (comMPI is not None)): #OC16042021 (i.e. several independent MPI calculations required)
        #if((_n_mpi > 1) and (_com_mpi is None) and (comMPI is not None)): #OC01032021 (i.e. several independent MPI calculations required)
            if((_opt_bl is None) or (_det is not None)): #OC01032021 (i.e. final mesh is known)
                nxny = 0
                if(_det is not None):
                    nxny = _det.nx*_det.ny
                else:
                    nxny = _mesh.nx*_mesh.ny

                itStart = 0; itEnd = nxny - 1 #? #OC03102021
                itStartEnd = [itStart, itEnd] #OC03102021 #Start and End indexes for aligned conjugated coordinate

            #OC03102021 (moved the lines below to be executed for _n_mpi > 1)
            fpCore = 'part_mi'
            fpExt = '.dat'
            if(_file_path is not None):
                iDot = _file_path.rfind('.')
                lenFP = len(_file_path)
                if((iDot >= 0) and (iDot >= lenFP - 5)):
                    fpCore = _file_path[0:iDot]
                    fpExt = _file_path[iDot:lenFP]

            nProcPerCalc = int(nProc/_n_mpi)
            nProcExtra = nProc - nProcPerCalc*_n_mpi

            arProcInGroup = [nProcPerCalc]*_n_mpi
            for i in range(nProcExtra): arProcInGroup[i] += 1
            
            #DEBUG
            #print('Rank #', rank,' arProcInGroup=', arProcInGroup)
            #sys.stdout.flush()
            #END DEBUG
                
            arRankMasters = [0]*_n_mpi #OC30102021
            rankMaster = -1

            rankStart = 0 #Start and end rank numbers in current Group
            rankEnd = -1 #0 #OC22042021
            iiGr = -1 #OC30102021
            for iGr in range(_n_mpi):
                arRankMasters[iGr] = rankStart #OC30102021
                rankEnd += arProcInGroup[iGr] #- 1 #OC22042021

                if((rank <= rankEnd) and (rankMaster < 0)): #OC30102021
                    rankMaster = rankStart #Master rank in current Group
                    iiGr = iGr #OC30102021

                #if(rank <= rankEnd): break
                rankStart = rankEnd + 1

            if(iiGr >= 0): iGr = iiGr

            #rankMaster = rankStart #OC30102021 (moved up, into the loop) #Master rank in current Group
            if(rankMaster < 0): rankMaster = 0 #OC30102021

                #OC25042021: The lines below were programmed for CSD calculaiton by parts (in attempt to save memory).
                # nfHalfCSD = 0.5*nxny*(nxny + 1)
                # nfNodesCSDperCalc = nfHalfCSD/_n_mpi
                # kfCSD = 0
                # itStart = 0
                # for ip in range(iGr + 1):
                #     kfCSD += nfNodesCSDperCalc
                #     itEnd = int(ceil(0.5*(sqrt(9 + 8*kfCSD) - 3.)))
                #     it_mi_1 = itEnd - 1
                #     kAux = 0.5*it_mi_1*(it_mi_1 + 3)
                #     iAux = kfCSD - kAux - 1 #OC16042021
                #     #iAux = kCSD - kAux - 1
                #     rat = iAux/(itEnd + 1)
                #     if(rat < 0.5): itEnd -= 1
                #     if(itEnd < itStart): itEnd = itStart
                #     #itStartEnd = [itStart, itEnd]
                #     if(ip < iGr): itStart = itEnd + 1

                #OC25042021: The method we currently use does full CSD calculation by each "MPI group":
                #itStart = 0; itEnd = nxny - 1 #? #OC03102021 (moved up)

            #DEBUG
            #print('Rank #', rank,' arProcInGroup=', arProcInGroup)
            #sys.stdout.flush()
            #END DEBUG

            nProc = arProcInGroup[iGr] #OC28102021 (moved from below)

            #OC25042021: Taking into account averaging
            #nPartTotOrig = _n_part_tot #OC28102021
            _n_part_tot = (int)(_n_part_tot/_n_mpi) #OC02112021: this is important here if not doPropCM

            # arNumModesForGroups = [_n_part_tot]*_n_mpi

            # nPartTotExtra = nPartTotOrig - _n_part_tot*_n_mpi #OC28102021
            # for iii in range(nPartTotExtra):
            #     arNumModesForGroups[iii] += 1
            #     arNumModesForGroups[nPartTotExtra - iii - 1] -= 1

            # for jjj in range(_n_mpi):
            #     jGrToDecrModes = _n_mpi - jjj - 1
            #     nModesForCurGroup = arNumModesForGroups[jGrToDecrModes]
            #     nModePerProcForCurGroup = int(float(nModesForCurGroup)/(arProcInGroup[jGrToDecrModes] - 1))
            #     nModeInExcessForCurGroup = nModesForCurGroup - nModePerProcForCurGroup*(arProcInGroup[jGrToDecrModes] - 1)

            #     if(nModeInExcessForCurGroup > 0):
            #         cntModeInExcessForCurGroup = 0
            #         for iGrToAddModes in range(jGrToDecrModes):
            #             nTotInGrToAddModes = arNumModesForGroups[iGrToAddModes]
            #             nProcInGrToAddModes = arProcInGroup[iGrToAddModes]
            #             nModePerProcForGrToAddModes = int(float(nTotInGrToAddModes)/(nProcInGrToAddModes - 1))
            #             nModeCanBeAdded = nTotInGrToAddModes - (nProcInGrToAddModes - 1)*nModePerProcForGrToAddModes

            #             if(nModeCanBeAdded <= 0): nModeCanBeAdded = nProcInGrToAddModes - 1
            #             allModesAdded = False
            #             for iModesAdded in range(nModeCanBeAdded):
            #                 arNumModesForGroups[iGrToAddModes] += 1
            #                 arNumModesForGroups[jGrToDecrModes] -= 1

            #                 cntModeInExcessForCurGroup += 1
            #                 if(cntModeInExcessForCurGroup == nModeInExcessForCurGroup): 
            #                     allModesAdded = True
            #                     break
            #             if(allModesAdded): break

            # _n_part_tot = arNumModesForGroups[iGr] #Number of Modes to be calculated by this MPI group

            # nModesPerWorkProc = int(_n_part_tot/(nProc - 1))
            # nExtraModesToRedist = _n_part_tot - nModesPerWorkProc*(nProc - 1)
            # if(nExtraModesToRedist > 0):
            #     nPossibleExtraModesInGr = nProc - 1 - nExtraModesToRedist

            # if(nPartTotExtra > 0): #OC28102021

            #     nMPItoAddModes = int(float(nPartTotExtra)/(nProc - 1) + 1.e-10)
            #     iMPI = int(float(rankMaster)/nProc + 1.e-10)
            #     if(iMPI < nMPItoAddModes): 
            #         _n_part_tot += nProc - 1
            #         #if(rank == rankMaster): _n_part_tot += nProc - 1
            #         #elif((rankMaster < rank) and (rank < (rankMaster + nProc))): _n_part_tot += 1
            #         #if(rank == nMPItoAddModes): 
            #         #_n_part_tot += nMPItoAddModes
            #     elif(iMPI == nMPItoAddModes):
            #         nRemainPartTotExtra = nPartTotExtra - nMPItoAddModes*(nProc - 1)
            #         if(nRemainPartTotExtra > 0): #_n_part_tot += 1
            #             _n_part_tot += nRemainPartTotExtra
            #             #if(rank == rankMaster): _n_part_tot += nRemainPartTotExtra
            #             #elif(rank - rankMaster <= nRemainPartTotExtra): _n_part_tot += 1
            # else:
            #     nTestPartPerProc = int(float(_n_part_tot)/(nProc - 1))
            #     nPartToRedistrFromOneMPI = _n_part_tot - nTestPartPerProc*(nProc - 1)

            #     if((nPartToRedistrFromOneMPI > 0) and (_n_mpi > 1)):
            #         #arExtraPartInMPI = [0]*_n_mpi
            #         #maxExtraPartCanBeAdded = nProc - 1 - 
            #         #for iMPI in range(nPartToRedistrFromOneMPI):

            #         if(float(rankMaster)/nProc < nPartToRedistrFromOneMPI): _n_part_tot += 1
            #         else: _n_part_tot -= 1

                #if(iMPI < nPartTotExtra): _n_part_tot += 1
                
                #for iMPI in range(_n_mpi): #OC28102021
                #    if(iMPI < nPartTotExtra): 
                #        if(rank - rankMaster < nProc): _n_part_tot += 1
                
                #itStartEnd = [itStart, itEnd] #OC03102021 (moved up) #Start and End indexes for aligned conjugated coordinate

            #nProc = arProcInGroup[iGr] #OC28102021 (moved up)
            
            #DEBUG
            #if(rank == rankMaster):
            #    print('Rank #', rank,' is Master; there are', nProc, 'processes in its group')
            #    sys.stdout.flush()
            #END DEBUG

            fp_cm = srwl_wfr_fn(_file_path, _type=7, _form='hdf5') #OC28062021
            _file_path = fpCore + '_' + repr(iGr) + fpExt #OC28062021: to change this (the input file name should not be changed!)

    else: #if not ((_char == 6) or (_char == 61) or (_char == 7)): #OC28042022
        if(_n_mpi > 1): raise Exception("Calculation with more than one \"group\" of MPI processes is is not supported for this radiation characteristic.")

    #DEBUG
    #if(rank == rankMaster):
    #    print('Rank #', rank,' is Master; there are', nProc, 'processes in its group')
    #    sys.stdout.flush()
    #END DEBUG

    #if(nProc <= 1): #OC050214
    #    _n_part_avg_proc = _n_part_tot

    iModeStart = 0 if(_ms <= 0) else int(_ms) #OC221122

    #OC04112020
    doPropCM = False
    if isinstance(_mag, list):
        if(len(_mag) > 0):
            if isinstance(_mag[0], SRWLWfr):
                doPropCM = True
                _mesh = copy(_mag[iModeStart].mesh) #OC23112022
                #_mesh = copy(_mag[0].mesh) #OC07112020

                #OC22092022
                if(_e_beam is not None):
                    m1 = _e_beam.partStatMom1
                    m1w = _mag[iModeStart].partBeam.partStatMom1 #OC23112022
                    #m1w = _mag[0].partBeam.partStatMom1
                    dz = m1.z - m1w.z
                    if(dz != 0.):
                        ebmNew = deepcopy(_e_beam)
                        ebmNew.drift(-dz) #propagating e-beam to the longitudinal position specified in the list of CMs 
                        _e_beam = ebmNew

    #OC30052017 (commented-out)
    #wfr = SRWLWfr() #Wavefronts to be used in each process
    #wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    #wfr.mesh.set_from_other(_mesh)    
    #wfr.partBeam = deepcopy(_e_beam)

    #arPrecParSR = [_sr_meth, _sr_rel_prec, 0, 0, 50000, 0, _sr_samp_fact] #to add npTraj, useTermin ([4], [5]) terms as input parameters
    arPrecParSR = [_sr_meth, _sr_rel_prec, 0, 0, 50000, 1, _sr_samp_fact] #to add npTraj, useTermin ([4], [5]) terms as input parameters

    #meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) if not doPropCM else None #OC04112020
    meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) #OC15112020 #OC30052017 (uncommented) #to ensure correct final mesh if _opt_bl is None
    #meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) if(_det is None) else _det.get_mesh() #OC06122016
    #Note: the case ((_det is not None) and (doMutual > 0)) is not supported

    file_path1 = copy(_file_path) #OC30052017
    file_path2 = None
    file_path_deg_coh1 = None #07052018
    file_path_deg_coh2 = None
    file_pathA = None
    meshRes2 = None #OC30052017
    meshResA = None #OC23122018
    
    #if(_char == 4): #Mutual Intensity Cuts vs X & Y are required
    if((_char == 4) or (_char == 40)): #OC02052018 #Mutual Intensity Cuts vs X & Y are required
        #meshRes2 = copy(meshRes)
        if(_char == 4): meshRes2 = copy(meshRes) #OC02052018
        
        #file_path1 += '.1'
        #file_path2 = copy(_file_path) + '.2'
        #file_path_deg_coh1 = copy(_file_path) + '.dc.1'
        #file_path_deg_coh2 = copy(_file_path) + '.dc.2'
        #OC19122018
        file_path1 = srwl_wfr_fn(_file_path, 31) # += '.1'
        file_path2 = srwl_wfr_fn(_file_path, 32) #copy(_file_path) + '.2'
        file_path_deg_coh1 = srwl_wfr_fn(_file_path, 41) #copy(_file_path) + '.dc.1'
        file_path_deg_coh2 = srwl_wfr_fn(_file_path, 42) #copy(_file_path) + '.dc.2'

    if((_char == 5) or (_char == 41)): #OC15072019 #Degree of Coherence Cuts vs X & Y are required
        if(_char == 5): meshRes2 = copy(meshRes) #OC15072019
        file_path_deg_coh1 = srwl_wfr_fn(_file_path, 41) #copy(_file_path) + '.dc.1'
        file_path_deg_coh2 = srwl_wfr_fn(_file_path, 42) #copy(_file_path) + '.dc.2'

    if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021
    #if(_char == 6): #OC02022021
        file_path1 = srwl_wfr_fn(_file_path, 3)
        #print('file_path1=', file_path1) #DEBUG

    if(_pres_ang == 2): #OC23122018 #Both coordinate and angular presentation characteristics are necessary
        #if((_char == 0) or (_char == 1) or (_char == 40)):
        if((_char == 0) or (_char == 1) or (_char == 40) or (_char == 41)): #OC15072019
            file_pathA = srwl_wfr_fn(_file_path, 2)

    wfr = SRWLWfr() #Wavefronts to be used in each process
    wfr2 = None #OC30052017
    wfrA = None #OC18062021 #OC23122018 (Average Wavefront data, e.g. for CSD / CMD calculations)
    
    if((_opt_bl is None) and (doMutual > 0)): #OC30052017
        if(_char == 2): #Cut vs X
            meshRes.ny = 1
            meshRes.yStart = _y0
            meshRes.yFin = _y0
        elif(_char == 3): #Cut vs Y
            meshRes.nx = 1
            meshRes.xStart = _x0
            meshRes.xFin = _x0
        #elif(_char == 4): #Cuts of Mutual Intensity vs X & Y
        elif((_char == 4) or (_char == 5)): #15072019 #Cuts of Mutual Intensity vs X & Y
            meshRes.ny = 1
            meshRes.yStart = _y0
            meshRes.yFin = _y0
            meshRes2.nx = 1
            meshRes2.xStart = _x0
            meshRes2.xFin = _x0
            wfr2 = SRWLWfr()
        #if((_char == 40) and (_pres_ang == 2)): #Mutual Intensity, Degree of Coherence, and Intensity if the Coordinate and Angular representation

    #OC04112020
    if not doPropCM:
        #OC30052017
        wfr.allocate(meshRes.ne, meshRes.nx, meshRes.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh.set_from_other(meshRes)    
        wfr.partBeam = deepcopy(_e_beam)
        if(wfr2 is not None):
            wfr2.allocate(meshRes2.ne, meshRes2.nx, meshRes2.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
            wfr2.mesh.set_from_other(meshRes2)    
            wfr2.partBeam = deepcopy(_e_beam)

    if(_det is not None): #OC06122016
        eStartOld = meshRes.eStart
        eFinOld = meshRes.eFin
        neOld = meshRes.ne

        meshRes = _det.get_mesh() 

        #OC17102021 (added)
        if(meshRes.eStart <= 0.): meshRes.eStart = eStartOld
        if(meshRes.eFin <= 0.): meshRes.eFin = eFinOld
        if(meshRes.ne < 1): meshRes.ne = neOld

    ePhIntegMult = 1
    if(_e_ph_integ == 1): #Integrate over Photon Energy
        eAvg = 0.5*(_mesh.eStart + _mesh.eFin)
        ePhIntegMult = 1000*(_mesh.eFin - _mesh.eStart)/eAvg #To obtain photon energy integrated Intensity in [ph/s/mm^2] assuming monochromatic Spectral Intensity in [ph/s/.1%bw/mm^2]
        wfr.mesh.eStart = eAvg 
        wfr.mesh.eFin = eAvg
        wfr.mesh.ne = 1
        if(wfr2 is not None):
            wfr2.mesh.eStart = eAvg #OC30052017
            wfr2.mesh.eFin = eAvg
            wfr2.mesh.ne = 1
        meshRes.eStart = eAvg #to check compatibility with _mesh_res is not None
        meshRes.eFin = eAvg
        meshRes.ne = 1
        if(meshRes2 is not None): #OC30052017
            meshRes2.eStart = eAvg #to check compatibility with _mesh_res is not None
            meshRes2.eFin = eAvg
            meshRes2.ne = 1

    calcSpecFluxSrc = False
    if(((_char == 10) or (_char == 11)) and (_mesh.nx == 1) and (_mesh.ny == 1)): #OC16042020
    #if((_char == 10) and (_mesh.nx == 1) and (_mesh.ny == 1)):
        calcSpecFluxSrc = True
        ePhIntegMult *= 1.e+06*(_mesh.xFin - _mesh.xStart)*(_mesh.yFin - _mesh.yStart) #to obtain Flux from Intensity (Flux/mm^2)

        #DEBUG
        #print('Const. for Flux computation was calculated: ePhIntegMult=', ePhIntegMult)
        #END DEBUG

    if not doPropCM: #OC04112020
        elecX0 = _e_beam.partStatMom1.x
        elecXp0 = _e_beam.partStatMom1.xp
        elecY0 = _e_beam.partStatMom1.y
        elecYp0 = _e_beam.partStatMom1.yp
        elecGamma0 = _e_beam.partStatMom1.gamma
        elecE0 = elecGamma0*(0.51099890221e-03) #Assuming electrons 
    
        elecSigXe2 = _e_beam.arStatMom2[0] #<(x-x0)^2>
        elecMXXp = _e_beam.arStatMom2[1] #<(x-x0)*(xp-xp0)>
        elecSigXpe2 = _e_beam.arStatMom2[2] #<(xp-xp0)^2>
        elecSigYe2 =_e_beam.arStatMom2[3] #<(y-y0)^2>
        elecMYYp = _e_beam.arStatMom2[4] #<(y-y0)*(yp-yp0)>
        elecSigYpe2 = _e_beam.arStatMom2[5] #<(yp-yp0)^2>
        elecRelEnSpr = sqrt(_e_beam.arStatMom2[10]) #<(E-E0)^2>/E0^2
        elecAbsEnSpr = elecE0*elecRelEnSpr
        #print('DEBUG MESSAGE: elecAbsEnSpr=', elecAbsEnSpr)
        #Consider taking into account other 2nd order moments?
    
        multX = 0.5/(elecSigXe2*elecSigXpe2 - elecMXXp*elecMXXp)
        BX = elecSigXe2*multX
        GX = elecSigXpe2*multX
        AX = elecMXXp*multX
        SigPX = 1/sqrt(2*GX)
        SigQX = sqrt(GX/(2*(BX*GX - AX*AX)))
        multY = 0.5/(elecSigYe2*elecSigYpe2 - elecMYYp*elecMYYp)
        BY = elecSigYe2*multY
        GY = elecSigYpe2*multY
        AY = elecMYYp*multY
        SigPY = 1/sqrt(2*GY)
        SigQY = sqrt(GY/(2*(BY*GY - AY*AY)))

        if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021
        #if(_char == 6): #OC18062021
            wfrA = SRWLWfr() #Average Wavefront params are required for eventual CSD / CMD calculation
            wfrA.xc = elecX0
            wfrA.yc = elecY0

    #_sr_rel_prec = int(_sr_rel_prec)
    
    _n_part_tot = int(_n_part_tot)
    _n_part_avg_proc = int(_n_part_avg_proc)
    if(_n_part_avg_proc <= 0): _n_part_avg_proc = 1
    _n_save_per = int(_n_save_per)

    nPartPerProc = _n_part_tot
    nSentPerProc = 0
    actNumPartTot = _n_part_tot #OC31102021
    #actNumPartTot = None #OC21012021
    
    if(nProc <= 1):
        if doPropCM: #OC26042022
            nModes = len(_mag)
            if(_n_part_tot > nModes): _n_part_tot = nModes
            nPartPerProc = _n_part_tot
            actNumPartTot = _n_part_tot
            #_n_part_avg_proc = 1

        _n_part_avg_proc = 1 #OC26042022 (??)
        #_n_part_avg_proc = _n_part_tot

    else: #OC050214: adjustment of all numbers of points, to make sure that sending and receiving are consistent

        if(not doPropCM):
             
            nPartPerProc = int(round(_n_part_tot/(nProc - 1)))
            nSentPerProc = int(round(nPartPerProc/_n_part_avg_proc)) #Number of sending acts made by each worker process
            if(nSentPerProc <= 0): #OC160116
                nSentPerProc = 1
                _n_part_avg_proc = nPartPerProc

            nPartPerProc = _n_part_avg_proc*nSentPerProc #Number of electrons treated by each worker process

            #DEBUG
            #print('Rank #', rank,': nProc=', nProc, 'nPartPerProc=', nPartPerProc, 'nSentPerProc=', nSentPerProc, '_n_part_avg_proc=', _n_part_avg_proc, 'iGr=', iGr)
            #sys.stdout.flush()
            #END DEBUG

        else: #OC19112020

            #OC30102021
            nModesTotToCalc = len(_mag)
            if(nModesTotToCalc > nPartTotOrig): nModesTotToCalc = nPartTotOrig

            arModesToCalcByThisWorker = []
            arModesToCalcByGroups = [0]*_n_mpi

            rankCnt = 0
            iCurGr = 0
            for im in range(nModesTotToCalc):
                if(arRankMasters is not None):
                    if(rankCnt in arRankMasters): 
                        iCurGr = arRankMasters.index(rankCnt)
                        rankCnt += 1
                        
                if(rankCnt == rank): 
                    arModesToCalcByThisWorker.append(im + iModeStart) #OC23112022
                    #arModesToCalcByThisWorker.append(im)
                
                arModesToCalcByGroups[iCurGr] += 1
                rankCnt += 1
                if(rankCnt >= nProcTot): rankCnt = 0

            _n_part_tot = arModesToCalcByGroups[iGr]

            nPartPerProc = len(arModesToCalcByThisWorker)
            nPartTot = _n_part_tot if(_n_part_tot < nModesTotToCalc) else nModesTotToCalc #OC30102021

            # nModes = len(_mag) #OC08082021
            # nPartTot = _n_part_tot if(_n_part_tot < nModes) else nModes #OC08082021
            
            # nPartPerProc = int(nPartTot/(nProc - 1)) #OC08082021
            # nRemain = nPartTot - nPartPerProc*(nProc - 1) #OC08082021

            # #nPartPerProc = int(_n_part_tot/(nProc - 1))
            # #nRemain = _n_part_tot - nPartPerProc*(nProc - 1)
            
            # if(nRemain > 0):
            #     if(rank - rankMaster <= nRemain): nPartPerProc += 1 #OC02032021
            #     #if(rank <= nRemain): nPartPerProc += 1

            nSentPerProc = int(nPartPerProc/_n_part_avg_proc)
            actNumPartTot = nPartTot #OC02032021
            #actNumPartTot = _n_part_tot
            
            if(nSentPerProc <= 0): #OC160116
                if(rank != rankMaster): #OC31102021
                    nSentPerProc = 1
                    _n_part_avg_proc = nPartPerProc
                else:
                    nPartPerProcEst = int(round(_n_part_tot/(nProc - 1)))
                    nSentPerProc = int(round(nPartPerProcEst/_n_part_avg_proc)) #Number of sending acts made by each worker process

            #DEBUG
            #print('Rank #', rank,': arModesToCalcByThisWorker=', arModesToCalcByThisWorker)
            #sys.stdout.flush()
            #END DEBUG
         
    #DEBUG
    #print('Rank #', rank,': nProc=', nProc, 'nPartPerProc=', nPartPerProc, 'nSentPerProc=', nSentPerProc, '_n_part_avg_proc=', _n_part_avg_proc, 'actNumPartTot=', actNumPartTot, 'iGr=', iGr)
    #sys.stdout.flush()
    #END DEBUG

    if(rank == rankMaster):
        log_path = srwl_uti_save_stat_wfr_emit_prop_multi_e_init(nProc, _n_part_avg_proc, actNumPartTot, iGr) #OC31102021

    #log_path, total_num_of_particles = srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, nProc, nPartPerProc, nSentPerProc, _n_part_avg_proc, actNumPartTot, iGr) #OC02032021
    #log_path, total_num_of_particles = srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, nProc, nPartPerProc, nSentPerProc, _n_part_avg_proc, actNumPartTot) #OC21012021
    #log_path, total_num_of_particles = srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, nProc, nPartPerProc, nSentPerProc, _n_part_avg_proc) #MR20012017

    #DEBUG
    #print('Rank #', rank,': _n_part_tot=', _n_part_tot, 'nPartPerProc=', nPartPerProc, 'nSentPerProc=', nSentPerProc, '_n_part_avg_proc=', _n_part_avg_proc)
    #sys.stdout.flush()
    #END DEBUG

    useGsnBmSrc = False
    usePtSrc = False #OC16102017
    if(isinstance(_mag, SRWLGsnBm)):
        useGsnBmSrc = True
        arPrecParSR = [_sr_samp_fact]
        _mag = deepcopy(_mag)
        _mag.x = elecX0
        _mag.xp = elecXp0
        _mag.y = elecY0
        _mag.yp = elecYp0
        #print('Gaussian Beam')
        #sys.exit()
    elif(isinstance(_mag, SRWLPtSrc)): 
        usePtSrc = True #OC16102017
        arPrecParSR = [_sr_samp_fact]
        _mag = deepcopy(_mag)
        _mag.x = elecX0
        _mag.xp = 0
        _mag.y = elecY0
        _mag.yp = 0

    resStokes = None
    workStokes = None
    #workStokes1cmp = None #OC04032024 (commented-out) #OC29122023
    resStokes2 = None #OC30052017
    workStokes2 = None
    #arAuxResSt12 = None #OC31052017
    #arAuxWorkSt12 = None #OC31052017
    arAuxResSt = None #OC24122018
    arAuxWorkSt = None #OC24122018

    resStokes3 = None #OC03052018
    workStokes3 = None
    #arAuxResSt123 = None #OC03052018
    #arAuxWorkSt123 = None #OC03052018

    resStokesA = None #OC23122018
    workStokesA = None #OC23122018

    iAvgProc = 0
    iSave = 0

    #doMutual = 0
    #if((_char >= 2) and (_char <= 4)): doMutual = 1

    #depTypeME_Approx = 3 #vs x&y (horizontal and vertical positions or angles); to be used only at _me_approx > 0
    #OC15092017
    depTypeInt = 3 #vs x&y (horizontal and vertical positions or angles); to be used only at _me_approx > 0
    #phEnME_Approx = _mesh.eStart
    
    #phEnInt = _mesh.eStart if not doPropCM else _mag[0].mesh.eStart #OC04112020
    phEnInt = _mesh.eStart #OC15092017
    #if(_me_approx == 1): #OC15092017
    #if((_mesh.nx <= 1) or (_mesh.ny <= 1)):
    if((_me_approx == 1) and ((_mesh.nx <= 1) or (_mesh.ny <= 1))): #OC01102017
        raise Exception("This calculation method requires more than one observation point in the horizontal and vertical directions")

    if(_mesh.ne > 1):
        #OC16042020
        if((_mesh.nx <= 1) and (_mesh.ny <= 1)): depTypeInt = 0 #vs e (photon energy or time)
        else:
            #OC15092017
            depTypeInt = 6 #vs e&x&y (photon energy or time, horizontal and vertical positions or angles)

        phEnInt = 0.5*(_mesh.eStart + _mesh.eFin)
        #depTypeME_Approx = 6 #vs e&x&y (photon energy or time, horizontal and vertical positions or angles)
        #phEnME_Approx = 0.5*(_mesh.eStart + _mesh.eFin)

    resEntityName = 'Intensity' #OC26042016
    resEntityNameDC = 'Degree of Coherence' #OC12072019
    
    resEntityUnits = 'ph/s/.1%bw/mm^2'
    if(_e_ph_integ > 0): resEntityUnits = 'ph/s/mm^2' #OC12072019
    resEntityUnitsDC = '' #OC12072019
    
    resLabelsToSave = ['Photon Energy', 'Horizontal Position', 'Vertical Position', resEntityName] #OC26042016
    resUnitsToSave = ['eV', 'm', 'm', resEntityUnits] #OC26042016

    resEntityNameA = None #OC24122018
    resEntityUnitsA = None
    resLabelsToSaveA = None
    resUnitsToSaveA = None
    
    if(_pres_ang == 1): #OC20112017
        resEntityName = 'Ang. Intensity Distr.'
        resEntityUnits = 'ph/s/.1%bw/mrad^2'
        if(_e_ph_integ > 0): resEntityUnits = 'ph/s/mrad^2' #OC12072019
        resLabelsToSave = ['Photon Energy', 'Horizontal Angle', 'Vertical Angle', resEntityName]
        resUnitsToSave = ['eV', 'rad', 'rad', resEntityUnits]
    elif(_pres_ang == 2): #OC24122018
        resEntityNameA = 'Ang. Intensity Distr.'
        resEntityUnitsA = 'ph/s/.1%bw/mrad^2'
        if(_e_ph_integ > 0): resEntityUnitsA = 'ph/s/mrad^2' #OC12072019
        resLabelsToSaveA = ['Photon Energy', 'Horizontal Angle', 'Vertical Angle', resEntityName]
        resUnitsToSaveA = ['eV', 'rad', 'rad', resEntityUnits]

    if(calcSpecFluxSrc == True):
        resEntityName = 'Flux'
        resEntityUnits = 'ph/s/.1%bw'
        if(_e_ph_integ > 0): resEntityUnits = 'ph/s' #OC12072019
        resLabelsToSave[3] = resEntityName #OC27032018
        resUnitsToSave[3] = resEntityUnits #OC27032018

    resLabelsToSaveDC = None #OC12072019
    resUnitsToSaveDC = None #OC12072019
    if((_char == 4) or (_char == 5) or (_char == 40) or (_char == 41)): #OC12072019 (Degree of Coherence required)
        resLabelsToSaveDC = copy(resLabelsToSave)
        resLabelsToSaveDC[3] = "Degree of Coherence"
        resUnitsToSaveDC = copy(resUnitsToSave)
        resUnitsToSaveDC[3] = ''

    #OC13112020 (moved this up)
    numComp = 1
    if((_char == 1) or (_char == 11) or (_char == 20)): numComp = 4

    if(_opt_bl is None): arPrecParSR[6] = 0 #Ensure non-automatic choice of numbers of points if there is no beamline

    bkpFileToBeSaved = False #OC14082018
    RxAvg = 0; RyAvg = 0; xcAvg = 0; ycAvg = 0 #OC21052020

    #resLabelsToSaveMutualHorCut = [resLabelsToSave[0], resLabelsToSave[1], 'Conj. ' + resLabelsToSave[1], 'Mutual ' + resLabelsToSave[3]] #OC03052018
    #resLabelsToSaveMutualVerCut = [resLabelsToSave[0], resLabelsToSave[2], 'Conj. ' + resLabelsToSave[2], 'Mutual ' + resLabelsToSave[3]]

    #if(((rank == 0) or (nProc == 1)) and (_opt_bl != None)): #calculate once the central wavefront in the master process (this has to be done only if propagation is required)
    #if(((rank == 0) or (nProc == 1)) and (_det is None)): #12/01/2017
    #OC23122018: need to better undestand the above change
    if(((rank == rankMaster) or (nProc == 1)) and (_opt_bl is not None) and (_det is None)): #OC02032021 
    #if(((rank == 0) or (nProc == 1)) and (_opt_bl is not None) and (_det is None)): #OC26022021 
    
    #Calculation of propagation of "central" beam / mode is required for obtaining mesh, but only if it is not known (for all workers)
    #OC26022021 Note: if _opt_bl is None then the final resulting mesh is known without any preliminary calculation

        if(useGsnBmSrc):
            srwl.CalcElecFieldGaussian(wfr, _mag, arPrecParSR)
            if(wfr2 is not None): srwl.CalcElecFieldGaussian(wfr2, _mag, arPrecParSR) #OC30052017
            #print('DEBUG: Commented-out: CalcElecFieldGaussian')
        elif(usePtSrc): #OC16102017
            srwl.CalcElecFieldPointSrc(wfr, _mag, arPrecParSR)
            if(wfr2 is not None): srwl.CalcElecFieldPointSrc(wfr2, _mag, arPrecParSR)

        elif(doPropCM): #OC04112020
            wfr = _mag[iModeStart] #OC221122
            #wfr = _mag[0] #OC13112020
            #wfr = deepcopy(_mag[0]) #OC12112020

            #OC10092022
            if(_e_beam is not None):
                m1 = _e_beam.partStatMom1; m1w = wfr.partBeam.partStatMom1
                dx = m1.x - m1w.x; dxp = m1.xp - m1w.xp; dy = m1.y - m1w.y; dyp = m1.yp - m1w.yp
                if((dx != 0.) or (dxp != 0.) or (dy != 0.) or (dyp != 0.)): 
                    wfr.sim_src_offset(_dx = dx, _dxp = dxp, _dy = dy, _dyp = dyp, _move_mesh=False, _copy=False)
            #OC01082022
            #if((wfr.partBeam.partStatMom1.x != _e_beam.partStatMom1.x) or (wfr.partBeam.partStatMom1.xp != _e_beam.partStatMom1.xp) or 
            #   (wfr.partBeam.partStatMom1.y != _e_beam.partStatMom1.y) or (wfr.partBeam.partStatMom1.yp != _e_beam.partStatMom1.yp)):
            #    wfr.sim_src_offset(_dx = (_e_beam.partStatMom1.x - wfr.partBeam.partStatMom1.x), _dxp = (_e_beam.partStatMom1.xp - wfr.partBeam.partStatMom1.xp), 
            #                       _dy = (_e_beam.partStatMom1.y - wfr.partBeam.partStatMom1.y), _dyp = (_e_beam.partStatMom1.yp - wfr.partBeam.partStatMom1.yp), _move_mesh=False, _copy=False)
                #wfr = wfr.sim_src_offset(_dx = (_e_beam.partStatMom1.x - wfr.partBeam.partStatMom1.x), _dxp = (_e_beam.partStatMom1.xp - wfr.partBeam.partStatMom1.xp), 
                #                         _dy = (_e_beam.partStatMom1.y - wfr.partBeam.partStatMom1.y), _dyp = (_e_beam.partStatMom1.yp - wfr.partBeam.partStatMom1.yp), _move_mesh=False, _copy=True)
            #DEBUG
            #wfr = _mag[4]
            #END DEBUG
            
        else:

            #print('Single-electron SR calculation ... ', end='') #DEBUG
            #t0 = time.time(); #DEBUG

            #print('DEBUG: wfr.mesh:', wfr.mesh)
            srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR)

            if(wfr2 is not None):
                #print('DEBUG: wfr2.mesh:', wfr2.mesh)
                srwl.CalcElecFieldSR(wfr2, 0, _mag, arPrecParSR) #OC30052017
            
            #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
            #print('DEBUG MESSAGE: CalcElecFieldSR called (rank:', rank,')')
            
        #print('DEBUG MESSAGE: Central Wavefront calculated')

        #print('Wavefront propagation calculation ... ', end='') #DEBUG
        #t0 = time.time(); #DEBUG

        #if(_w_wr != 0.): #OC26032016
        if(_wr != 0.): #OC07092016
            wfr.Rx = _wr
            wfr.Ry = _wr
            if(wfr2 is not None): #OC30052017
                wfr2.Rx = _wr
                wfr2.Ry = _wr

        if(_wre > 0.): #OC05012017
            wfr.dRx = _wre
            wfr.dRy = _wre
            if(wfr2 is not None): #OC30052017
                wfr2.dRx = _wre
                wfr2.dRy = _wre
        
        if(_opt_bl is not None): #OC16012017
            #if(_rand_opt): _opt_bl.randomize() #OC24042020: don't randomize here yet

            #DEBUG
            #print('Rank #', rank, 'Before Central Propag.:')
            #print('wfr.mesh.nx=', wfr.mesh.nx, 'wfr.mesh.ny=', wfr.mesh.ny)
            #print('wfr.mesh.xStart=', wfr.mesh.xStart, 'wfr.mesh.xFin=', wfr.mesh.xFin, 'wfr.mesh.yStart=', wfr.mesh.yStart, 'wfr.mesh.yFin=', wfr.mesh.yFin)
            #sys.stdout.flush()
            #END DEBUG

            srwl.PropagElecField(wfr, _opt_bl)

        #DEBUG
        #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
        ##meshRes.set_from_other(wfr.mesh) #DEBUG
        #meshRes = copy(wfr.mesh) #DEBUG
        #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual) #DEBUG
        #print('Rank #', rank, 'After Central Propag.: nx=', meshRes.nx, ' ny=', meshRes.ny) #DEBUG
        #print('meshRes.xStart=', meshRes.xStart, 'meshRes.xFin=', meshRes.xFin, 'meshRes.yStart=', meshRes.yStart, 'meshRes.yFin=', wfr.mesh.yFin)
        #sys.stdout.flush()
        #wfr.calc_stokes(resStokes) #DEBUG
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual) #DEBUG
        #sys.exit() #DEBUG
        #print('wavefront / mode propagated') #DEBUG
        #END DEBUG

        #print('DEBUG: Commented-out: PropagElecField')
        #print('DEBUG MESSAGE: Central Wavefront propagated')
        #if(_pres_ang > 0):
        if(_pres_ang == 1): #OC23122018
            
            wfr.unitElFldAng = 1 #OC20112017 (to ensure result in [ph/s/.1%bw/mrad^2])
            srwl.SetRepresElecField(wfr, 'a')
            if(wfr2 is not None):
                wfr2.unitElFldAng = 1 #OC20112017
                srwl.SetRepresElecField(wfr2, 'a') #OC30052017
            #print('DEBUG: Commented-out: SetRepresElecField')

        #meshRes.set_from_other(wfr.mesh)
        #if(_det is None): meshRes.set_from_other(wfr.mesh) #OC06122016
        #else: meshRes = _det.get_mesh() #??

        #OC04112020
        if not doPropCM: meshRes.set_from_other(wfr.mesh)
        else: meshRes = copy(wfr.mesh)
        #meshRes.set_from_other(wfr.mesh) #OC12012016

        if(_pres_ang == 2): #OC23122018
            wfr.unitElFldAng = 1 #?
            srwl.SetRepresElecField(wfr, 'a')

            if(meshResA is None): meshResA = SRWLRadMesh()
            meshResA.set_from_other(wfr.mesh)

        #if(wfr2 is not None): meshRes2.set_from_other(wfr2.mesh) #OC30052017

        if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021
        #if(_char == 6): #OC18062021
            if(wfrA is None): wfrA = SRWLWfr()
            wfrA.numTypeElFld = wfr.numTypeElFld
            wfrA.Rx = wfr.Rx
            wfrA.Ry = wfr.Ry
            wfrA.dRx = wfr.dRx
            wfrA.dRy = wfr.dRy
            wfrA.xc = elecX0
            wfrA.yc = elecY0
            wfrA.avgPhotEn = wfr.avgPhotEn
            wfrA.presCA = wfr.presCA
            wfrA.presFT = wfr.presFT
            wfrA.unitElFld = wfr.unitElFld
            wfrA.unitElFldAng = wfr.unitElFldAng
            wfrA.mesh = wfr.mesh #OC28062021

        if(doMutual > 0):
            if(_char == 2): #Cut vs X
                meshRes.ny = 1
                meshRes.yStart = _y0
                meshRes.yFin = _y0
            elif(_char == 3): #Cut vs Y
                meshRes.nx = 1
                meshRes.xStart = _x0
                meshRes.xFin = _x0
            #elif(_char == 4): #Cuts of Mutual Intensity vs X & Y
            elif((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                meshRes.ny = 1
                meshRes.yStart = _y0
                meshRes.yFin = _y0
                if(wfr2 is None): meshRes2.set_from_other(wfr.mesh) #OC31052017
                else: meshRes2.set_from_other(wfr2.mesh) #OC31052017
                meshRes2.nx = 1
                meshRes2.xStart = _x0
                meshRes2.xFin = _x0

        #OC13112020 (moved this up)
        if(_char == 41): #If Intensity and Degree of Coherence is required, send over average Wavefront Radii of Curvature for eventual subtraction of Quadratic Phase Terms (to improve accuracy of Degrre of Coherence)
            RxAvg = wfr.Rx
            if(0.2*abs(RxAvg) < abs(wfr.dRx)): RxAvg = 0
            RyAvg = wfr.Ry
            if(0.2*abs(RyAvg) < abs(wfr.dRy)): RyAvg = 0
            #arMesh.extend((RxAvg, RyAvg, wfr.xc, wfr.yc))

        if(nProc > 1): #send resulting mesh to all workers
            #comMPI.send(wfr.mesh, dest=)
            arMesh = array('f', [meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny])
            
            #arMesh2 = None #OC30052017
            #if(_char == 4): #Cuts of Mutual Intensity vs X & Y
            if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                arMesh2 = array('f', [meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny])
                for ii in range(len(arMesh2)): #OC31052017
                    arMesh.append(arMesh2[ii])

            if(meshResA is not None): #OC23122018
                arMeshA = array('f', [meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny])
                for ii in range(len(arMeshA)): arMesh.append(arMeshA[ii])

            #comMPI.Bcast([arMesh, MPI.FLOAT], root=MPI.ROOT)
            #comMPI.Bcast([arMesh, MPI.FLOAT])

            #OC13112020 (moved this up)
            #OC21052020
            if(_char == 41): #If Intensity and Degree of Coherence is required, send over average Wavefront Radii of Curvature for eventual subtraction of Quadratic Phase Terms (to improve accuracy of Degrre of Coherence)
                #RxAvg = wfr.Rx
                #if(0.2*abs(RxAvg) < abs(wfr.dRx)): RxAvg = 0
                #RyAvg = wfr.Ry
                #if(0.2*abs(RyAvg) < abs(wfr.dRy)): RyAvg = 0
                arMesh.extend((RxAvg, RyAvg, wfr.xc, wfr.yc))

            #print('DEBUG MESSAGE: Master is about to broadcast mesh of Propagated central wavefront')
            for iRank in range(nProc - 1):
                dst = iRank + 1 + rankMaster #OC02032021
                #dst = iRank + 1
                #print("msg %d: sending data from %d to %d" % (iRank, rank, dst)) #an he

                comMPI.Send([arMesh, MPI.FLOAT], dest=dst)

                #DEBUG
                #print('Mesh data sent from Master to Worker #', dst)
                #END DEBUG

                #OC30052017
                #if(_char == 4): #Cuts of Mutual Intensity vs X & Y
                #    comMPI.Send([arMesh2, MPI.FLOAT], dest=dst)

            #DEBUG
            #print('Mesh of propagated central wavefront broadcasted')
            #sys.stdout.flush()
            #END DEBUG

        #DEBUG
        #print('meshRes: ne=', meshRes.ne, 'eStart=', meshRes.eStart, 'eFin=', meshRes.eFin)
        #print('Rank #', rank, 'doMutual=', doMutual, 'meshRes: nx=', meshRes.nx, 'xStart=', meshRes.xStart, 'xFin=', meshRes.xFin, 'ny=', meshRes.ny, 'yStart=', meshRes.yStart, 'yFin=', meshRes.yFin)
        #sys.stdout.flush()
        #sys.exit(0)
        #END DEBUG

        resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, _n_comp = numComp) #OC18072021
        #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

        if(doPropCM):
            wfr.calc_stokes(resStokes, _n_stokes_comp = numComp) #OC18072021
            #wfr.calc_stokes(resStokes) #OC13112020

            #DEBUG
            #print('Mode #0 propagated by Master; resStokes defined')
            #sys.stdout.flush()
            #END DEBUG

        #wfr.calc_stokes(resStokes) #OC190414: don't take into account the first "central" beam!
        #workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
        #OC06042017 (commented-out workStokes = ...)

        lenArSt0 = len(resStokes.arS) #OC24122018
        lenArSt = lenArSt0
        
        #if(_char == 4): #OC31052017 #Cuts of Mutual Intensity vs X & Y
        if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity or Degree of Coherence vs X & Y
            resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual, _n_comp = numComp) #OC18072021
            #resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
            #lenArSt12 = len(resStokes.arS) + len(resStokes2.arS)
            #arAuxResSt12 = array('f', [0]*lenArSt12)
            lenArSt += len(resStokes2.arS) #OC24122018

        #if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
        if((_char == 40) or (_char == 41)): #OC15072019 #Intensity and Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
            resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1, _n_comp = numComp) #OC18072021
            #resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)

            #DEBUG
            #print('Allocating resStokes3:', meshRes.yStart, meshRes.yFin, meshRes.ny)
            resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1, _n_comp = numComp) #OC18072021
            #resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)

            if(doPropCM): #OC13112020
                wfr.calc_stokes(resStokes2, _n_stokes_comp=numComp, _rx_avg=RxAvg, _ry_avg=RyAvg, _xc_avg=xcAvg, _yc_avg=ycAvg)
                wfr.calc_stokes(resStokes3, _n_stokes_comp=numComp, _rx_avg=RxAvg, _ry_avg=RyAvg, _xc_avg=xcAvg, _yc_avg=ycAvg)

                #DEBUG
                #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, _file_path + '.test_dcx0.dat', 1,
                #                           _arLabels = ['Photon Energy', '(x1+x2)/2','(x1-x2)/2','Degree of Coherence'], _arUnits = ['eV','m','m',''], _mutual=2, _cmplx=0)
                #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, _file_path + '.test_dcy0.dat', 1,
                #                           _arLabels = ['Photon Energy', '(y1+y2)/2','(y1-y2)/2','Degree of Coherence'], _arUnits = ['eV','m','m',''], _mutual=2, _cmplx=0)
                #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, _file_path + '.test_mix0.dat', 1,
                #                           _arLabels = ['Photon Energy', 'x1','x2','Mutual Intensity'], _arUnits = ['eV','m','m',''], _mutual=1, _cmplx=1)
                #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, _file_path + '.test_miy0.dat', 1,
                #                           _arLabels = ['Photon Energy', 'y1','y2','Mutual Intensity'], _arUnits = ['eV','m','m',''], _mutual=1, _cmplx=1)
                #sys.exit()
                #END DEBUG

            #lenArSt123 = len(resStokes.arS) + len(resStokes2.arS) + len(resStokes3.arS)
            #arAuxResSt123 = array('f', [0]*lenArSt123)
            lenArSt += (len(resStokes2.arS) + len(resStokes3.arS)) #OC24122018

        if(meshResA is not None): #OC23122018
            resStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny, _n_comp = numComp) #OC18072021
            #resStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny)
            lenArSt += len(resStokesA.arS)

        if((lenArSt > lenArSt0) and ((arAuxResSt is None) or (len(arAuxResSt) < lenArSt))): arAuxResSt = array('f', [0]*lenArSt) #OC24122018

        #iAvgProc += 1 #OC190414 (commented-out)
        #iSave += 1

    #slaves = [] #an he
    #print('DEBUG MESSAGE: rank=', rank, ' rankMaster=', rankMaster)

    #OC13112020 (moved this up)
    #numComp = 1
    #if((_char == 1) or (_char == 11) or (_char == 20)): numComp = 4 #OC16042020
    #if(_char == 20): numComp = 4 #OC16012017

    #if(_opt_bl is None): arPrecParSR[6] = 0 #Ensure non-automatic choice of numbers of points if there is no beamline
    #bkpFileToBeSaved = False #OC14082018
    #RxAvg = 0; RyAvg = 0; xcAvg = 0; ycAvg = 0 #OC21052020

    if((rank > rankMaster) or (nProc == 1)): #OC02032021
    #if((rank > 0) or (nProc == 1)):
        #if((nProc > 1) and (_opt_bl != None)): #receive mesh for the resulting wavefront from the master
        #if((nProc > 1) and (_opt_bl is not None) and (_det is None)): #OC12012017 #receive mesh for the resulting wavefront from the master
        #if((nProc > 1) and (_det is None)): #OC16012017 #receive mesh for the resulting wavefront from the master
        if((nProc > 1) and (_det is None) and (_opt_bl is not None)): #OC28022021 #receive mesh for the resulting wavefront from the master
           #In this part, mesh information is received by workers from master
            #arMesh = array('f', [0]*9)
            nNumToRecv = 9
            #if(_char == 4): nNumToRecv = 18 #OC31052017 #Cuts of Mutual Intensity vs X & Y
            if((_char == 4) or (_char == 5)): nNumToRecv = 18 #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
            if(_pres_ang == 2): nNumToRecv += 9 #OC26122018

            #If Intensity and Degree of Coherence is required, send over average Wavefront Radii of Curvature for eventual subtraction of Quadratic Phase Terms (to improve accuracy of Degrre of Coherence)
            if(_char == 41): nNumToRecv += 4 #OC21052020 
            
            arMesh = array('f', [0]*nNumToRecv)

            #_stat = MPI.Status() #an he
            #comMPI.Recv([arMesh, MPI.FLOAT], source=0)
            comMPI.Recv([arMesh, MPI.FLOAT], source=MPI.ANY_SOURCE)
            #comMPI.Bcast([arMesh, MPI.FLOAT], root=0)
            #print("received mesh %d -> %d" % (_stat.Get_source(), rank))
            
            meshRes.eStart = arMesh[0]
            meshRes.eFin = arMesh[1]
            meshRes.ne = int(arMesh[2])
            meshRes.xStart = arMesh[3]
            meshRes.xFin = arMesh[4]
            meshRes.nx = int(arMesh[5])
            meshRes.yStart = arMesh[6]
            meshRes.yFin = arMesh[7]
            meshRes.ny = int(arMesh[8])

            iStA = 9 #OC23122018

            #arMesh2 = None #OC30052017
            #if(_char == 4): #Cuts of Mutual Intensity vs X & Y
            if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                #arMesh2 = array('f', [0]*9)
                #_stat = MPI.Status() #an he
                #comMPI.Recv([arMesh2, MPI.FLOAT], source=MPI.ANY_SOURCE)
                #comMPI.Bcast([arMesh, MPI.FLOAT], root=0)
                #print("received mesh %d -> %d" % (_stat.Get_source(), rank))
                meshRes2.eStart = arMesh[9] #OC31052017
                meshRes2.eFin = arMesh[10]
                meshRes2.ne = int(arMesh[11])
                meshRes2.xStart = arMesh[12]
                meshRes2.xFin = arMesh[13]
                meshRes2.nx = int(arMesh[14])
                meshRes2.yStart = arMesh[15]
                meshRes2.yFin = arMesh[16]
                meshRes2.ny = int(arMesh[17])
                iStA += 9 #OC21052020
                #iStA = 18 #OC23122018

            if(_pres_ang == 2): #OC23122018
                if(meshResA is None):
                    meshResA = SRWLRadMesh()
                
                meshResA.eStart = arMesh[iStA]
                meshResA.eFin = arMesh[iStA + 1]
                meshResA.ne = int(arMesh[iStA + 2])
                meshResA.xStart = arMesh[iStA + 3]
                meshResA.xFin = arMesh[iStA + 4]
                meshResA.nx = int(arMesh[iStA + 5])
                meshResA.yStart = arMesh[iStA + 6]
                meshResA.yFin = arMesh[iStA + 7]
                meshResA.ny = int(arMesh[iStA + 8])
                iStA += 9 #OC21052020 

            if(_char == 41): #OC21052020
                RxAvg = arMesh[iStA]
                RyAvg = arMesh[iStA + 1]
                xcAvg = arMesh[iStA + 2]
                ycAvg = arMesh[iStA + 3]

            #DEBUG
            #print('Mesh data from Master received by Worker #', rank)
            #sys.stdout.flush()
            #END DEBUG

        #sys.exit(0) #DEBUG

        #OC17022021 (commented-out the following, because nRadPt*, nStPt* doesn't seem to be used)
        #nRadPt = meshRes.ne*meshRes.nx*meshRes.ny
        #if(doMutual > 0): nRadPt *= nRadPt
        #nStPt = nRadPt*4

        #nRadPt2 = None #OC30052017
        #nStPt2 = None
        ##if(_char == 4): #Cuts of Mutual Intensity vs X & Y
        #if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
        #    nRadPt2 = meshRes2.ne*meshRes2.nx*meshRes2.ny
        #    if(doMutual > 0): nRadPt2 *= nRadPt2 #OC03052018
        #    nStPt2 = nRadPt2*4

        #nRadPt3 = None #OC03052018
        #nStPt3 = None
        ##if(_char == 40): #Intensity and Cuts of Mutual Intensity vs X & Y
        #if((_char == 40) or (_char == 41)): #OC15072019 #Intensity and Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
        #    nRadPt2 = meshRes.ne*meshRes.nx
        #    nRadPt2 *= nRadPt2
        #    nStPt2 = nRadPt2*4
        #    nRadPt3 = meshRes.ne*meshRes.ny
        #    nRadPt3 *= nRadPt3
        #    nStPt3 = nRadPt3*4

        #DEBUG
        #print('Rank #', rank, ': meshRes.nx=', meshRes.nx, 'meshRes.ny=', meshRes.ny)
        #print('Rank #', rank, ': nRadPt=', nRadPt, 'nRadPt2=', nRadPt2, 'nRadPt3=', nRadPt3)
        #print('Rank #', rank, ': nStPt=', nStPt, 'nStPt2=', nStPt2, 'nStPt3=', nStPt3)
        #print('Rank #', rank, ': nPartPerProc=', nPartPerProc)
        #sys.stdout.flush()
        #END DEBUG
        
        randAr = array('d', [0]*6) #for random Gaussian numbers

        #random.seed(rank) #old

        #OCTEST OC10062021
        #random.seed(4*123) #Deterministic seeding
        #print('Artificial test seed (as for rank #4)')
        #END OCTEST

        random.seed(rank*123) #Deterministic seeding

        #DEBUG OC16102021
        #random.seed((rank + 1230)*123) #Deterministic seeding for debug
        #random.seed((rank + 336)*123) #Deterministic seeding for debug
        #END DEBUG

        #random.seed() #Pseudo-random seeding based on instant time

        #OC24042021 (commented-out the lines below)
        #if(_n_mpi > 1): #OC02032021 (to ensure independent calculations in each MPI group)
        #    random.seed()

        newSeed = random.randint(0, 1000000)

        #DEBUG
        #print('Rank #', rank, ': newSeed=', newSeed)
        #sys.stdout.flush()
        #END DEBUG

        random.seed(newSeed)

        iAuxSendCount = 0 #for DEBUG

        #OC18022021
        arElFldToSend = None #To be used only in the case on nProc > 1 and _char == 6
        lenHalfArToSend = 0; lenArToSend = 0
        #To allocate array for sending Electric Field, using the numbers of points in the final mesh (to be taken from _det or from first test propagated wavefront)
        if(((_char == 6) or (_char == 61) or (_char == 7)) and (nProc > 1)): #OC20062021
        #if((_char == 6) and (nProc > 1)): #OC27022021
        #if(_char == 6): 
            #Asumes that meshRes is already defined at this point for workers
            lenHalfArToSend = meshRes.ne*meshRes.nx*meshRes.ny*2 #for arEx and arEy (consider introducing some logic of arEx or arEy is not used)
            lenArToSend = 2*lenHalfArToSend
            arElFldToSend = array('f', [0]*lenArToSend)
        
        for i in range(nPartPerProc): #loop over macro-electrons

            if((_me_approx == 0) and (not doPropCM)): #OC04112020
            #if(_me_approx == 0): #OC05042017 #General method
                
                if(_rand_meth == 1):
                    for ir in range(5): #to expend to 6D eventually
                        randAr[ir] = random.gauss(0, 1)
                elif(_rand_meth == 2):
                    if(nProc > 1):
                        #iArg = i*(nProc - 1) + rank - rankMaster #OC02032021: not necessary, since rank is unique
                        iArg = i*(nProc - 1) + rank
                        a1 = srwl_uti_math_seq_halton(iArg, 2)
                        a2 = srwl_uti_math_seq_halton(iArg, 3)
                        a3 = srwl_uti_math_seq_halton(iArg, 5)
                        a4 = srwl_uti_math_seq_halton(iArg, 7)
                        a5 = srwl_uti_math_seq_halton(iArg, 11) #?
                    elif(nProc == 1):
                        i_p_1 = i + 1
                        a1 = srwl_uti_math_seq_halton(i_p_1, 2)
                        a2 = srwl_uti_math_seq_halton(i_p_1, 3)
                        a3 = srwl_uti_math_seq_halton(i_p_1, 5)
                        a4 = srwl_uti_math_seq_halton(i_p_1, 7)
                        a5 = srwl_uti_math_seq_halton(i_p_1, 11) #?
                    twoPi = 2*pi
                    twoPi_a2 = twoPi*a2
                    twoPi_a4 = twoPi*a4
                    m2_log_a1 = -2.0*log(a1)
                    m2_log_a3 = -2.0*log(a3)
                    randAr[0] = sqrt(m2_log_a1)*cos(twoPi_a2)
                    randAr[1] = sqrt(m2_log_a1)*sin(twoPi_a2)
                    randAr[2] = sqrt(m2_log_a3)*cos(twoPi_a4)
                    randAr[3] = sqrt(m2_log_a3)*sin(twoPi_a4)
                    randAr[4] = sqrt(m2_log_a1)*cos(twoPi*a3) #or just random.gauss(0,1) depends on cases #why not using a5?
                    randAr[5] = a5
                elif(_rand_meth == 3):
                    #to program LPtau sequences here
                    continue

                #DEBUG
                #if(i == 0):
                #    randAr = array('d', [0,0,0,0,0])
                #if(i == 1):
                #    randAr = array('d', [0,0,0,-2,0])
                #END DEBUG

                auxPXp = SigQX*randAr[0]
                auxPX = SigPX*randAr[1] + AX*auxPXp/GX
                wfr.partBeam.partStatMom1.x = elecX0 + auxPX
                wfr.partBeam.partStatMom1.xp = elecXp0 + auxPXp
                auxPYp = SigQY*randAr[2]
                auxPY = SigPY*randAr[3] + AY*auxPYp/GY
                wfr.partBeam.partStatMom1.y = elecY0 + auxPY
                wfr.partBeam.partStatMom1.yp = elecYp0 + auxPYp
                #wfr.partBeam.partStatMom1.gamma = (elecEn0 + elecAbsEnSpr*randAr[4])/0.51099890221e-03 #Relative Energy
                #wfr.partBeam.partStatMom1.gamma = elecGamma0*(1 + elecAbsEnSpr*randAr[4]/elecE0)
                wfr.partBeam.partStatMom1.gamma = elecGamma0*(1 + elecRelEnSpr*randAr[4]) #OC28122016

                if(wfr2 is not None): #OC30052017
                    wfr2.partBeam.partStatMom1.x = elecX0 + auxPX
                    wfr2.partBeam.partStatMom1.xp = elecXp0 + auxPXp
                    wfr2.partBeam.partStatMom1.y = elecY0 + auxPY
                    wfr2.partBeam.partStatMom1.yp = elecYp0 + auxPYp
                    wfr2.partBeam.partStatMom1.gamma = elecGamma0*(1 + elecRelEnSpr*randAr[4]) #OC28122016

                #Consider taking into account other 2nd order moments?

            elif((_me_approx == 1) and (not doPropCM)): #OC04112020
            #elif(_me_approx == 1): #OC05042017 #Numerical integration only over electron energy

                if(_rand_meth == 1):
                    randAr[0] = random.gauss(0, 1)
                elif(_rand_meth == 2):
                    if(nProc > 1):
                        iArg = i*(nProc - 1) + rank
                        a1 = srwl_uti_math_seq_halton(iArg, 2)
                        a2 = srwl_uti_math_seq_halton(iArg, 3)
                        #a3 = srwl_uti_math_seq_halton(iArg, 5)
                        #a4 = srwl_uti_math_seq_halton(iArg, 7)
                        #a5 = srwl_uti_math_seq_halton(iArg, 11) #?
                    elif(nProc == 1):
                        i_p_1 = i + 1
                        a1 = srwl_uti_math_seq_halton(i_p_1, 2)
                        a2 = srwl_uti_math_seq_halton(i_p_1, 3)
                        #a3 = srwl_uti_math_seq_halton(i_p_1, 5)
                        #a4 = srwl_uti_math_seq_halton(i_p_1, 7)
                        #a5 = srwl_uti_math_seq_halton(i_p_1, 11) #?
                    twoPi = 2*pi
                    twoPi_a2 = twoPi*a2
                    #twoPi_a4 = twoPi*a4
                    m2_log_a1 = -2.0*log(a1)
                    #m2_log_a3 = -2.0*log(a3)
                    randAr[0] = sqrt(m2_log_a1)*cos(twoPi_a2)
                    #randAr[1] = sqrt(m2_log_a1)*sin(twoPi_a2)
                    #randAr[2] = sqrt(m2_log_a3)*cos(twoPi_a4)
                    #randAr[3] = sqrt(m2_log_a3)*sin(twoPi_a4)
                    #randAr[4] = sqrt(m2_log_a1)*cos(twoPi*a3) #or just random.gauss(0,1) depends on cases #why not using a5?
                    #randAr[5] = a5
                elif(_rand_meth == 3):
                    #to program LPtau sequences here
                    continue

                wfr.partBeam.partStatMom1.x = elecX0
                wfr.partBeam.partStatMom1.xp = elecXp0
                wfr.partBeam.partStatMom1.y = elecY0
                wfr.partBeam.partStatMom1.yp = elecYp0
                wfr.partBeam.partStatMom1.gamma = elecGamma0*(1 + elecRelEnSpr*randAr[0]) #OC05042017
                if(wfr2 is not None): #OC30052017
                    wfr2.partBeam.partStatMom1.x = elecX0
                    wfr2.partBeam.partStatMom1.xp = elecXp0
                    wfr2.partBeam.partStatMom1.y = elecY0
                    wfr2.partBeam.partStatMom1.yp = elecYp0
                    wfr2.partBeam.partStatMom1.gamma = elecGamma0*(1 + elecRelEnSpr*randAr[0])

            if not doPropCM: #OC04112020
                #OC06042017 (added for _me_approx == 1 and possibly other future methods)
                wfr.partBeam.arStatMom2[0] = elecSigXe2 #<(x-x0)^2>
                wfr.partBeam.arStatMom2[1] = elecMXXp #<(x-x0)*(xp-xp0)>
                wfr.partBeam.arStatMom2[2] = elecSigXpe2 #<(xp-xp0)^2>
                wfr.partBeam.arStatMom2[3] = elecSigYe2 #<(y-y0)^2>
                wfr.partBeam.arStatMom2[4] = elecMYYp #<(y-y0)*(yp-yp0)>
                wfr.partBeam.arStatMom2[5] = elecSigYpe2 #<(yp-yp0)^2>
                wfr.partBeam.arStatMom2[10] = elecRelEnSpr*elecRelEnSpr #<(E-E0)^2>/E0^2
                if(wfr2 is not None): #OC30052017
                    wfr2.partBeam.arStatMom2[0] = elecSigXe2 #<(x-x0)^2>
                    wfr2.partBeam.arStatMom2[1] = elecMXXp #<(x-x0)*(xp-xp0)>
                    wfr2.partBeam.arStatMom2[2] = elecSigXpe2 #<(xp-xp0)^2>
                    wfr2.partBeam.arStatMom2[3] = elecSigYe2 #<(y-y0)^2>
                    wfr2.partBeam.arStatMom2[4] = elecMYYp #<(y-y0)*(yp-yp0)>
                    wfr2.partBeam.arStatMom2[5] = elecSigYpe2 #<(yp-yp0)^2>
                    wfr2.partBeam.arStatMom2[10] = elecRelEnSpr*elecRelEnSpr #<(E-E0)^2>/E0^2

            #Consider taking into account other 2nd order moments?

            #reset mesh, because it may be modified by CalcElecFieldSR and PropagElecField
            #print('Numbers of points (before re-setting): nx=', wfr.mesh.nx, ' ny=', wfr.mesh.ny) #DEBUG
            curWfrMesh = wfr.mesh #OC02042016
            newWfrMesh = _mesh if(_opt_bl is not None) else meshRes #OC16012017

            #if((curWfrMesh.ne != _mesh.ne) or (curWfrMesh.nx != _mesh.nx) or (curWfrMesh.ny != _mesh.ny)):
            #    wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny)
            if((curWfrMesh.ne != newWfrMesh.ne) or (curWfrMesh.nx != newWfrMesh.nx) or (curWfrMesh.ny != newWfrMesh.ny)):  #OC16012017
                #DEBUG
                #print('curWfrMesh: nx=', curWfrMesh.nx, ' ny=', curWfrMesh.ny)
                #print('newWfrMesh: nx=', newWfrMesh.nx, ' ny=', newWfrMesh.ny)
                #END DEBUG
                wfr.allocate(newWfrMesh.ne, newWfrMesh.nx, newWfrMesh.ny)

            if(_opt_bl is None): wfr.mesh.set_from_other(meshRes) #OC16012017
            else: wfr.mesh.set_from_other(_mesh)

            if(wfr2 is not None): #OC30052017
                curWfrMesh2 = wfr2.mesh
                newWfrMesh2 = _mesh if(_opt_bl is not None) else meshRes2

                if((curWfrMesh2.ne != newWfrMesh2.ne) or (curWfrMesh2.nx != newWfrMesh2.nx) or (curWfrMesh2.ny != newWfrMesh2.ny)):
                    wfr2.allocate(newWfrMesh2.ne, newWfrMesh2.nx, newWfrMesh2.ny)

                if(_opt_bl is None): wfr2.mesh.set_from_other(meshRes2)
                else: wfr2.mesh.set_from_other(_mesh)

            if(_e_ph_integ == 1):
                if(_rand_meth == 1):
                    ePh = random.uniform(_mesh.eStart, _mesh.eFin)
                else:
                    ePh = _mesh.eStart + (_mesh.eFin - _mesh.eStart)*randAr[5]
                    
                wfr.mesh.eStart = ePh
                wfr.mesh.eFin = ePh
                wfr.mesh.ne = 1
                if(wfr2 is not None):
                    wfr2.mesh.eStart = ePh
                    wfr2.mesh.eFin = ePh
                    wfr2.mesh.ne = 1

            wfr.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
            wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
            if(wfr2 is not None):
                wfr2.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
                wfr2.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time

            if(nProc == 1):
                if(useGsnBmSrc): print('i=', i, 'Gaussian Beam Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp)
                elif(usePtSrc): print('i=', i, 'Point Source Coord.: x=', wfr.partBeam.partStatMom1.x, 'y=', wfr.partBeam.partStatMom1.y)
                #elif(doPropCM): print('Mode:', iMode) #OC09112020
                #elif(doPropCM): print('Mode:', i) #OC04112020
                elif(doPropCM): print('Mode:', i + iModeStart) #OC221122
                else: print('i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
                
                if(_e_ph_integ == 1): print('Eph=', wfr.mesh.eStart)

            #DEBUG OC16102021
            #print('Rank:', rank, 'i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
            #print('DEBUG: re-defining macro-electron initial conditions: OLD values:')
            #print('i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
            #wfr.partBeam.partStatMom1.x = -3.662765020071964e-05 #OC: these initial conditions lead to wfr propagation off the final mesh in Example20 - to check the propagation issues(!)
            #wfr.partBeam.partStatMom1.xp = 4.2604539754324834e-05
            #wfr.partBeam.partStatMom1.y = -5.014957889048704e-06
            #wfr.partBeam.partStatMom1.yp = 2.4359620246785347e-06
            #wfr.partBeam.partStatMom1.gamma = 3.004478203792218/0.51099890221e-03
            #print('DEBUG: re-defining macro-electron initial conditions: NEW values:')
            #print('i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
            #sys.stdout.flush()
            #END DEBUG

            if(calcSpecFluxSrc): #consider taking into account _rand_meth != 1 here
                xObs = random.uniform(_mesh.xStart, _mesh.xFin)
                wfr.mesh.xStart = xObs
                wfr.mesh.xFin = xObs
                yObs = random.uniform(_mesh.yStart, _mesh.yFin)
                wfr.mesh.yStart = yObs
                wfr.mesh.yFin = yObs
                #DEBUG
                #print('xObs=', xObs, 'yObs=', yObs)
                #END DEBUG

            iMode = i + iModeStart #OC23112022
            #iMode = i #OC19112020
            try:
                if(useGsnBmSrc):
                    _mag.x = wfr.partBeam.partStatMom1.x
                    _mag.xp = wfr.partBeam.partStatMom1.xp
                    _mag.y = wfr.partBeam.partStatMom1.y
                    _mag.yp = wfr.partBeam.partStatMom1.yp
                    srwl.CalcElecFieldGaussian(wfr, _mag, arPrecParSR)
                    if(wfr2 is not None): srwl.CalcElecFieldGaussian(wfr2, _mag, arPrecParSR) #OC30052017
                    #print('DEBUG: Commented-out: CalcElecFieldGaussian')
                    #print('Gaussian wavefront calc. done')
                elif(usePtSrc):
                    _mag.x = wfr.partBeam.partStatMom1.x
                    _mag.xp = 0
                    _mag.y = wfr.partBeam.partStatMom1.y
                    _mag.yp = 0
                    srwl.CalcElecFieldPointSrc(wfr, _mag, arPrecParSR)
                    if(wfr2 is not None): srwl.CalcElecFieldPointSrc(wfr2, _mag, arPrecParSR) #OC30052017

                elif(doPropCM): #OC04112020

                    #wfr = _mag[i]

                    if(nProc > 1): iMode = arModesToCalcByThisWorker[i] #OC26042022
                    #iMode = arModesToCalcByThisWorker[i] #OC30102021
                    #if(nProc > 1): iMode = (nProc - 1)*(_n_mpi*i + int(float(rankMaster)/float(nProc) + 1.e-10)) + rank - rankMaster - 1 #OC27102021
                    #if(nProc > 1): iMode = (nProc - 1)*i + rank - 1 #OC19112020
                    #if(nProc > 1): iMode = (rank - 1)*nPartPerProc + i #OC19112020

                    #DEBUG
                    #print('rank:', rank, ' iMode=', iMode, 'assigned')
                    #sys.stdout.flush()
                    #END DEBUG
                    
                    wfr = _mag[iMode] #OC19112020

                    #OC10092022
                    if(_e_beam is not None):
                        m1 = _e_beam.partStatMom1; m1w = wfr.partBeam.partStatMom1
                        dx = m1.x - m1w.x; dxp = m1.xp - m1w.xp; dy = m1.y - m1w.y; dyp = m1.yp - m1w.yp
                        if((dx != 0.) or (dxp != 0.) or (dy != 0.) or (dyp != 0.)): 
                            wfr.sim_src_offset(_dx = dx, _dxp = dxp, _dy = dy, _dyp = dyp, _move_mesh=False, _copy=False)
                    #OC01082022
                    #if((wfr.partBeam.partStatMom1.x != _e_beam.partStatMom1.x) or (wfr.partBeam.partStatMom1.xp != _e_beam.partStatMom1.xp) or 
                    #   (wfr.partBeam.partStatMom1.y != _e_beam.partStatMom1.y) or (wfr.partBeam.partStatMom1.yp != _e_beam.partStatMom1.yp)):
                    #    wfr.sim_src_offset(_dx = (_e_beam.partStatMom1.x - wfr.partBeam.partStatMom1.x), _dxp = (_e_beam.partStatMom1.xp - wfr.partBeam.partStatMom1.xp), 
                    #                       _dy = (_e_beam.partStatMom1.y - wfr.partBeam.partStatMom1.y), _dyp = (_e_beam.partStatMom1.yp - wfr.partBeam.partStatMom1.yp), _move_mesh=False, _copy=False)
                        #wfr = wfr.sim_src_offset(_dx = (_e_beam.partStatMom1.x - wfr.partBeam.partStatMom1.x), _dxp = (_e_beam.partStatMom1.xp - wfr.partBeam.partStatMom1.xp), 
                        #                         _dy = (_e_beam.partStatMom1.y - wfr.partBeam.partStatMom1.y), _dyp = (_e_beam.partStatMom1.yp - wfr.partBeam.partStatMom1.yp), _move_mesh=False, _copy=True)
                    #DEBUG
                    #print('Rank #', rank, ': Mode #', iMode, 'assigned')
                    #sys.stdout.flush()
                    #END DEBUG
                    
                else:

                    #print('Single-electron SR calculation ... ', end='') #DEBUG
                    #t0 = time.time(); #DEBUG
                    #print('arPrecParSR[6]=', arPrecParSR[6])

                    #DEBUG
                    #print('Numbers of points (after re-setting):') #DEBUG
                    #print('ne=', wfr.mesh.ne, 'eStart=', wfr.mesh.eStart, 'eFin=', wfr.mesh.eFin)
                    #print('nx=', wfr.mesh.nx, 'xStart=', wfr.mesh.xStart, 'xFin=', wfr.mesh.xFin)
                    #print('ny=', wfr.mesh.ny, 'yStart=', wfr.mesh.yStart, 'yFin=', wfr.mesh.yFin)
                    #print('zStart=', wfr.mesh.zStart)
                    #END DEBUG
                    
                    srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR) #calculate Electric Field emitted by current electron
                    if(wfr2 is not None): srwl.CalcElecFieldSR(wfr2, 0, _mag, arPrecParSR) #OC30052017

                    #DEBUG OC10102021
                    #if(rank == 203): 
                    #    print('    rank=', rank, 'iMode=', iMode, ': E-field calculated')
                    #    sys.stdout.flush()
                    #END DEBUG

                    #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                    #print('DEBUG: Commented-out: CalcElecFieldSR')
                    #print('DEBUG MESSAGE: CalcElecFieldSR called (rank:', rank,')')

                if(_opt_bl is not None):

                    #print('Wavefront propagation calculation ... ', end='') #DEBUG
                    #t0 = time.time(); #DEBUG

                    #if(_w_wr != 0.): #OC26032016
                    if(_wr != 0.): #OC07092016
                        wfr.Rx = _wr
                        wfr.Ry = _wr

                    if(_wre > 0.): #OC05012017
                        wfr.dRx = _wre
                        wfr.dRy = _wre

                    if(_rand_opt): _opt_bl.randomize() #OC24042020
                    #DEBUG OC24042020
                    #print('TEST: Coordinates of optical element center:', _opt_bl.arOpt[13].x, _opt_bl.arOpt[13].y)
                    #print('TEST: Coordinates of optical element center:', 0.5*(_opt_bl.arOpt[15].mesh.yStart + _opt_bl.arOpt[15].mesh.yFin))
                    #print('TEST: Coordinates of optical element center:', 0.5*(_opt_bl.arOpt[2].mesh.yStart + _opt_bl.arOpt[2].mesh.yFin))
                    #END DEBUG
                    #print('Before srwl.PropagElecField(wfr, _opt_bl)') #DEBUG

                    #DEBUG
                    #print('wfr.mesh.xStart=', wfr.mesh.xStart, 'wfr.mesh.xFin=', wfr.mesh.xFin, 'wfr.mesh.yStart=', wfr.mesh.yStart, 'wfr.mesh.yFin=', wfr.mesh.yFin)
                    #END DEBUG

                    if(not (doPropCM and (iMode == iModeStart) and (_det is None))): #OC23112022
                    #if(not (doPropCM and (iMode == 0) and (_det is None))): #OC17022021
                    #if(not (doPropCM and (iMode == 0) and (_det == 0))): #OC03022021
                    #if(not (doPropCM and (iMode == 0))): #OC19112020
                    #if(not (doPropCM and (i == 0))): #OC13112020

                        #DEBUG OC10102021
                        #if(rank == 203): 
                        #    print('    rank=', rank, 'iMode=', iMode, ': starting propag. E-field')
                        #    sys.stdout.flush()
                        #END DEBUG

                        srwl.PropagElecField(wfr, _opt_bl) #propagate Electric Field emitted by the electron

                        if(_det is not None): srwl.ResizeElecFieldMesh(wfr, meshRes, [0, 1]) #OC01082022

                        #DEBUG OC10102021
                        #if(rank == 203): 
                        #    print('    rank=', rank, 'iMode=', iMode, ': E-field calculated')
                        #    sys.stdout.flush()
                        #END DEBUG

                        #DEBUG
                        #print('Rank #', rank, ': mode #', iMode, 'propagated')
                        #sys.stdout.flush()
                        #END DEBUG

                    #DEBUG
                    #print('wfr.mesh.xStart=', wfr.mesh.xStart, 'wfr.mesh.xFin=', wfr.mesh.xFin, 'wfr.mesh.yStart=', wfr.mesh.yStart, 'wfr.mesh.yFin=', wfr.mesh.yFin)
                    #END DEBUG

                    #print('srwl.PropagElecField(wfr, _opt_bl) OK') #DEBUG
                    #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                    #print('DEBUG: Commented-out: PropagElecField')

                    #DEBUG
                    #if(i == 48):
                    #    arI1 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
                    #    srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
                    #    srwl_uti_save_intens_ascii(arI1, wfr.mesh, os.path.join(os.getcwd(), 'data_CDI', 'debug_int_pr_se.dat'))
                    #    sys.exit()
                    #END DEBUG

                #if(_pres_ang > 0):
                if(_pres_ang == 1): #OC23122018
                    
                    wfr.unitElFldAng = 1 #OC20112017 (to have result in [ph/s/.1%bw/mrad^2] vs [rad])
                    srwl.SetRepresElecField(wfr, 'a')
                    if(wfr2 is not None):
                        wfr2.unitElFldAng = 1 #OC20112017
                        srwl.SetRepresElecField(wfr2, 'a')

                    #print('DEBUG: Commented-out: SetRepresElecField')

            except:

                #DEBUG
                #if doPropCM:
                #    print('Rank=', rank, ' Mode=', iMode, ' i=', i, ' nPartPerProc=', nPartPerProc)
                #    sys.stdout.flush()
                #END DEBUG

                traceback.print_exc()

            #if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021 (moved down)
            ##if(_char == 6): #OC18062021
            #    if(i == 0): #OC26102021
            #        if(wfrA is None): wfrA = SRWLWfr()

            #        wfrA.numTypeElFld = wfr.numTypeElFld
            #        wfrA.Rx = wfr.Rx
            #        wfrA.Ry = wfr.Ry
            #        wfrA.dRx = wfr.dRx
            #        wfrA.dRy = wfr.dRy

            #        if doPropCM: #OC25102021
            #            wfrA.xc = wfr.xc
            #            wfrA.yc = wfr.yc
            #        else: #OC25102021: make sure if this is correct (probably the "if doPropCM" treatment should be done for all cases?)
            #            wfrA.xc = elecX0
            #            wfrA.yc = elecY0

            #        wfrA.avgPhotEn = wfr.avgPhotEn
            #        wfrA.presCA = wfr.presCA
            #        wfrA.presFT = wfr.presFT
            #        wfrA.unitElFld = wfr.unitElFld
            #        wfrA.unitElFldAng = wfr.unitElFldAng
            #        wfrA.mesh = wfr.mesh #OC28062021

            meshWork = deepcopy(wfr.mesh)
            meshWork2 = None #OC30052017
            meshWorkA = None #OC24122018
            if(doMutual > 0):
                if(_char == 2):
                    meshWork.ny = 1
                    meshWork.yStart = _y0
                    meshWork.yFin = _y0
                elif(_char == 3):
                    meshWork.nx = 1
                    meshWork.xStart = _x0
                    meshWork.xFin = _x0
                #elif(_char == 4): #OC30052017 #Cuts of Mutual Intensity vs X & Y
                elif((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                    meshWork.ny = 1
                    meshWork.yStart = _y0
                    meshWork.yFin = _y0
                    meshWork2 = deepcopy(wfr.mesh)
                    meshWork2.nx = 1
                    meshWork2.xStart = _x0
                    meshWork2.xFin = _x0

            #DEBUG
            #print('meshWork.xStart=', meshWork.xStart, 'meshWork.xFin=', meshWork.xFin, 'meshWork.nx=', meshWork.nx, 'meshWork.yStart=', meshWork.yStart, 'meshWork.yFin=', meshWork.yFin, 'meshWork.ny=', meshWork.ny)
            #print('meshRes.xStart=', meshRes.xStart, 'meshRes.xFin=', meshRes.xFin, 'meshRes.nx=', meshRes.nx, 'meshRes.yStart=', meshRes.yStart, 'meshRes.yFin=', meshRes.yFin, 'meshRes.ny=', meshRes.ny)
            #END DEBUG
            
            #OC06042017 (commented-out the above, entered workStokes = ... below)
            if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021
            #if(_char != 6): #OC03022021
                workStokes = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, meshWork.yStart, meshWork.yFin, meshWork.ny, doMutual, _n_comp = numComp) #OC18072021
                #workStokes = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, meshWork.yStart, meshWork.yFin, meshWork.ny, doMutual)

            #if(_char == 4): #Cuts of Mutual Intensity vs X & Y
            if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                workStokes2 = SRWLStokes(1, 'f', meshWork2.eStart, meshWork2.eFin, meshWork2.ne, meshWork2.xStart, meshWork2.xFin, meshWork2.nx, meshWork2.yStart, meshWork2.yFin, meshWork2.ny, doMutual, _n_comp = numComp) #OC18072021
                #workStokes2 = SRWLStokes(1, 'f', meshWork2.eStart, meshWork2.eFin, meshWork2.ne, meshWork2.xStart, meshWork2.xFin, meshWork2.nx, meshWork2.yStart, meshWork2.yFin, meshWork2.ny, doMutual)
                     
            #if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
            if((_char == 40) or (_char == 41)): #OC15072019 #Intensity and Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                workStokes2 = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, _y0, _y0, 1, _mutual=1, _n_comp = numComp) #OC18072021
                #workStokes2 = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, _y0, _y0, 1, _mutual=1)
                workStokes3 = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, _x0, _x0, 1, meshWork.yStart, meshWork.yFin, meshWork.ny, _mutual=1, _n_comp = numComp) #OC18072021
                #workStokes3 = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, _x0, _x0, 1, meshWork.yStart, meshWork.yFin, meshWork.ny, _mutual=1)

            if((_char == 6) or (_char == 61) or (_char == 7)): #OC03102021
            #if((_det is not None) and ((_char == 6) or (_char == 61) or (_char == 7))): #OC20062021 (consider doing this for other cases!)
            #if((_det is not None) and (_char == 6)): #OC03022021 (consider doing this for other cases!)

                #DEBUG OC10102021
                #if(rank == 203): 
                #    print('    rank=', rank, 'iMode=', iMode, ': staring resizing E-field to required final mesh')
                #    sys.stdout.flush()
                #END DEBUG

                #OC: if wfr.mesh and meshRes have the same basic params, perhaps calling the following function is not required(?)
                srwl.ResizeElecFieldMesh(wfr, meshRes, [0, 1])

                #if(i == 0): #OC26102021 (moved from top)
                if(wfrA is None): wfrA = SRWLWfr()

                if(wfrA.avgPhotEn <= 0): #I.e. if Average Wavefront was not filled-out

                    wfrA.numTypeElFld = wfr.numTypeElFld
                    wfrA.Rx = wfr.Rx
                    wfrA.Ry = wfr.Ry
                    wfrA.dRx = wfr.dRx
                    wfrA.dRy = wfr.dRy

                    if doPropCM: #OC25102021
                        wfrA.xc = wfr.xc
                        wfrA.yc = wfr.yc
                    else: #OC25102021: make sure if this is correct (probably the "if doPropCM" treatment should be done for all cases?)
                        wfrA.xc = elecX0
                        wfrA.yc = elecY0

                    wfrA.avgPhotEn = wfr.avgPhotEn
                    wfrA.presCA = wfr.presCA
                    wfrA.presFT = wfr.presFT
                    wfrA.unitElFld = wfr.unitElFld
                    wfrA.unitElFldAng = wfr.unitElFldAng
                    wfrA.mesh = copy(wfr.mesh) #OC27102021
                    #wfrA.mesh = wfr.mesh #OC28062021

                #DEBUG OC17102021
                #import numpy as np
                #print('rank=', rank, ': Trying to test E-field after Resizing at iMode=', iMode)
                #np.asarray_chkfinite(wfr.arEx, dtype=float)
                #np.asarray_chkfinite(wfr.arEy, dtype=float)
                #print('rank=', rank, 'iMode=', iMode, ': No NaN or Inf found in E-field')
                #sys.stdout.flush()
                #    print('    rank=', rank, 'iMode=', iMode, ': E-field resized to the final mesh')
                #    sys.stdout.flush()
                #END DEBUG

                #DEBUG
                #print('rank=', rank, 'iMode=', iMode, ': resizing to Detector Mesh done')
                #print('wfr.mesh.xStart=', wfr.mesh.xStart, 'wfr.mesh.xFin=', wfr.mesh.xFin, 'wfr.mesh.nx=', wfr.mesh.nx)
                #print('wfr.mesh.yStart=', wfr.mesh.yStart, 'wfr.mesh.yFin=', wfr.mesh.yFin, 'wfr.mesh.ny=', wfr.mesh.ny)
                #sys.stdout.flush()
                #END DEBUG

            #OC04022021 (moved from below)
            if(resStokes is None):
                #nComp = 4
                #if((_char == 6) or (_char == 61) or (_char == 7)): nComp = 1 #OC20062021
                ##if(_char == 6): nComp = 1 #OC04022021 ??
                #OC18072021 (commented-out the above, using numComp instead)

                if(not (((_char == 6) or (_char == 61) or (_char == 7)) and (nProc > 1))): #OC20062021
                #if(not ((_char == 6) and (nProc > 1))): #OC18022021
                    resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, numComp, itStartEnd) #OC18072021
                    #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, nComp, itStartEnd) #OC03032021
                    #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, nComp)

            if(_me_approx == 0): #OC05042017 #General case of numerical integration over 5D phase space of electron beam
                
                if(_char == 20): 
                    wfr.copy_comp(workStokes) #OC15012017: copy electric field components to Stokes structure
                elif(_char == 0): #OC15092017

                    #DEBUG
                    #print('About to define workStokes')
                    #END DEBUG
                    if(not (doPropCM and (iMode == iModeStart) and (_det is None))): #OC23112022
                    #if(not (doPropCM and (iMode == 0) and (_det is None))): #OC01082022
                    #if(not (doPropCM and (iMode == 0))): #OC19112020
                    #if(not (doPropCM and (i == 0))): #OC13112020
                        srwl.CalcIntFromElecField(workStokes.arS, wfr, 6, 0, depTypeInt, phEnInt, 0., 0.)
                    
                elif((_char == 1) or (_char == 11)): #OC16042020
                    #meshWorkStokes = workStokes.mesh #DEBUG
                    #print('workStokes: ne=', meshWorkStokes.ne, 'nx=', meshWorkStokes.nx, 'ny=', meshWorkStokes.ny, 'len(arS)=', len(workStokes.arS))
                    srwl.CalcIntFromElecField(workStokes.arS, wfr, -5, 0, depTypeInt, phEnInt, 0., 0.) #All Stokes

                elif((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021
                #elif(_char == 6): #OC03022021

                    #DEBUG
                    #print('resStokes.mesh.xStart=', resStokes.mesh.xStart, ' resStokes.mesh.xFin=', resStokes.mesh.xFin, ' resStokes.mesh.nx=', resStokes.mesh.nx)
                    #print('resStokes.mesh.yStart=', resStokes.mesh.yStart, ' resStokes.mesh.yFin=', resStokes.mesh.yFin, ' resStokes.mesh.ny=', resStokes.mesh.ny)
                    #print('depTypeInt=', depTypeInt, ' phEnInt=', phEnInt)
                    #END DEBUG

                    if(nProc == 1): #OC18022021
                        #Extract Intensity from Electric Field only at sequential execution (otherwise Electric Field should be sent directly to master)
                        #Calculate single-e mutual intensity from electric fied and add it to resStokes.arS (s0)
                        #NOTE: Common Quadratic Phase Terms can be subtracted before calculating the mutual intensity
                        intSumType = 1 #calculation of Intensity with instant averaging
                        if(doPropCM): intSumType = 2 #adding of new Intensity value to previous one

                        #DEBUG
                        #t0 = time.time()
                        #END DEBUG

                        arMethPar = [0]*20 #OC03032021 (in principle, this is not required if(nProc == 1) ?)
                        arMethPar[0] = intSumType; arMethPar[1] = i
                        if((_n_mpi > 1) and (itStartEnd is not None)):
                            arMethPar[18] = itStartEnd[0]
                            arMethPar[19] = itStartEnd[1]

                        #if(_opt_bl is not None): srwl.ResizeElecFieldMesh(wfr, resStokes.mesh, [0, 1]) #OC03102021 (changed lines before) #OC01102021 (added)

                        srwl.CalcIntFromElecField(resStokes.arS, wfr, -1, 8, depTypeInt, phEnInt, 0., 0., arMethPar) #OC03032021: this call is supposed to update / extract one main Stokes component
                        #srwl.CalcIntFromElecField(resStokes.arS, wfr, -1, 8, depTypeInt, phEnInt, 0., 0., [intSumType, i]) #OC03032021: this call is supposed to update / extract one main Stokes component
                        #srwl.CalcIntFromElecField(resStokes.arS, wfr, 6, 8, depTypeInt, phEnInt, 0., 0., [intSumType, i]) #One main Stokes component

                        #DEBUG
                        #print(resStokes.arS[0], resStokes.arS[resStokes.mesh.nx*resStokes.mesh.ny*2 + 2], resStokes.arS[(resStokes.mesh.nx*resStokes.mesh.ny*2)*2 + 2*2])
                        #END DEBUG

                        #DEBUG/TEST of a function:
                        #srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [5, -1, wfrA.Rx, wfrA.Ry, wfrA.xc, wfrA.yc]) #Subtract common quadratic phase terms from CSD
                        #END DEBUG/TEST of a function

                        #DEBUG
                        #print('CSD Update lasted:', round(time.time() - t0, 6), 's')
                        #sys.stdout.flush()
                        #END DEBUG

                    #DEBUG
                    #print('MI updated')
                    #END DEBUG
                    
                else:
                    #DEBUG
                    #print('About to define workStokes')
                    #END DEBUG
                    if(not (doPropCM and (iMode == iModeStart))): #OC23112022
                    #if(not (doPropCM and (iMode == 0))): #OC19112020
                    #if(not (doPropCM and (i == 0))): #OC13112020
                    
                        wfr.calc_stokes(workStokes, _n_stokes_comp=numComp) #calculate Stokes parameters from Electric Field

                        if(workStokes2 is not None): #OC30052017
                            #OC21052020
                            if(wfr2 is None): wfr.calc_stokes(workStokes2, _n_stokes_comp=numComp, _rx_avg=RxAvg, _ry_avg=RyAvg, _xc_avg=xcAvg, _yc_avg=ycAvg)
                            else: wfr2.calc_stokes(workStokes2, _n_stokes_comp=numComp, _rx_avg=RxAvg, _ry_avg=RyAvg, _xc_avg=xcAvg, _yc_avg=ycAvg)
                            #if(wfr2 is None): wfr.calc_stokes(workStokes2, _n_stokes_comp=numComp)
                            #else: wfr2.calc_stokes(workStokes2, _n_stokes_comp=numComp)

                        if(workStokes3 is not None): #OC03052018
                            wfr.calc_stokes(workStokes3, _n_stokes_comp=numComp, _rx_avg=RxAvg, _ry_avg=RyAvg, _xc_avg=xcAvg, _yc_avg=ycAvg) #OC21052020
                            #wfr.calc_stokes(workStokes3, _n_stokes_comp=numComp)

                if(_pres_ang == 2): #23122018
                    wfr.unitElFldAng = 1 #?
                    srwl.SetRepresElecField(wfr, 'a')
                    meshWorkA = deepcopy(wfr.mesh)
                    workStokesA = SRWLStokes(1, 'f', meshWorkA.eStart, meshWorkA.eFin, meshWorkA.ne, meshWorkA.xStart, meshWorkA.xFin, meshWorkA.nx, meshWorkA.yStart, meshWorkA.yFin, meshWorkA.ny, _n_comp = numComp) #OC18072021
                    #workStokesA = SRWLStokes(1, 'f', meshWorkA.eStart, meshWorkA.eFin, meshWorkA.ne, meshWorkA.xStart, meshWorkA.xFin, meshWorkA.nx, meshWorkA.yStart, meshWorkA.yFin, meshWorkA.ny)
                    srwl.CalcIntFromElecField(workStokesA.arS, wfr, 6, 0, depTypeInt, phEnInt, 0., 0.)

                    #DEBUG
                    #srwl_uti_save_intens_ascii(workStokesA.arS, workStokesA.mesh, copy(_file_path) + '.debug', 1)
                    #END DEBUG

            elif(_me_approx == 1): #OC05042017 #Numerical integration only over electron energy, convolution over transverse phase space

                if(_char == 0): #Total intensity
                    #DEBUG
                    #print('DEBUG: 2nd order e-beam moments after eventual propagation:')
                    #print('sigX=', sqrt(wfr.partBeam.arStatMom2[0]))
                    #print('mXXp=', wfr.partBeam.arStatMom2[1])
                    #print('sigXp=', sqrt(wfr.partBeam.arStatMom2[2]))
                    #print('sigY=', sqrt(wfr.partBeam.arStatMom2[3]))
                    #print('mYYp=', wfr.partBeam.arStatMom2[4])
                    #print('sigYp=', sqrt(wfr.partBeam.arStatMom2[5]))
                    #print('relEnSpr=', sqrt(wfr.partBeam.arStatMom2[10]))
                    #print('wfr.Rx=', wfr.Rx, ' wfr.Ry=', wfr.Ry)
                    #print('wfr.arElecPropMatr=', wfr.arElecPropMatr)
                    #END DEBUG

                    #srwl.CalcIntFromElecField(workStokes.arS, wfr, 6, 1, depTypeME_Approx, phEnME_Approx, _x0, _y0)
                    srwl.CalcIntFromElecField(workStokes.arS, wfr, 6, 1, depTypeInt, phEnInt, _x0, _y0) #OC30052017

                    #DEBUG
                    #arTest = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny)
                    ##srwl.CalcIntFromElecField(arTest, wfr, 6, 1, depTypeME_Approx, phEnME_Approx, _x0, _y0)
                    #srwl.CalcIntFromElecField(arTest, wfr, 6, 0, depTypeME_Approx, phEnME_Approx, _x0, _y0)

                    #srwl_uti_save_intens_ascii(arTest, wfr.mesh, _file_path + '.debug', numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)

                    ##srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path + '.debug', numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)
                    #if(i == 0): raise Exception("DEBUG: STOP") #OC06042017
                    ##srwl.CalcIntFromElecField(workStokes.arS, wfr, 6, 0, depTypeME_Approx, phEnME_Approx, _x0, _y0)
                    #END DEBUG

                elif(_char == 1): #Four Stokes components
                    
                    #OC29122023
                    srwl.CalcIntFromElecField(workStokes.arS, wfr, -5, 1, depTypeInt, phEnInt, 0., 0.) #All Stokes (implement extraction of Stokes components, with convolution, in CalcIntFromElecField)

                #if(_pres_ang == 2): #23122018
                #    #To make convolution taking into account angular divergancies only!
                #                
                #    wfr.unitElFldAng = 1 #?
                #    srwl.SetRepresElecField(wfr, 'a')
                #    meshWorkA = deepcopy(wfr.mesh)
                #    workStokesA = SRWLStokes(1, 'f', meshWorkA.eStart, meshWorkA.eFin, meshWorkA.ne, meshWorkA.xStart, meshWorkA.xFin, meshWorkA.nx, meshWorkA.yStart, meshWorkA.yFin, meshWorkA.ny)
                #    #srwl.CalcIntFromElecField(workStokesA.arS, wfr, 6, 1, depTypeInt, phEnInt, 0, 0)

            #DEBUG
            #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, 1)
            #END DEBUG

            #OC04022021 (moved up)
            #if(resStokes is None):
            #    resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
                #DEBUG
                #print('resStokes #2: ne=', resStokes.mesh.ne, 'eStart=', resStokes.mesh.eStart, 'eFin=', resStokes.mesh.eFin)
                #END DEBUG

            lenArSt0 = 0
            if(resStokes is not None): #OC18022021
                lenArSt0 = len(resStokes.arS)
            lenArSt = lenArSt0

            #if(_char == 4): #OC31052017 #Cuts of Mutual Intensity vs X & Y
            if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                if(resStokes2 is None):
                    resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual, _n_comp = numComp) #OC18072021
                    #resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
                #if(arAuxResSt12 is None):
                #    lenArSt12 = len(resStokes.arS) + len(resStokes2.arS)
                #    arAuxResSt12 = array('f', [0]*lenArSt12)
                lenArSt += len(resStokes2.arS) #OC24122018

            #if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
            if((_char == 40) or (_char == 41)): #OC15072019 #Intensity and Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
                if(resStokes2 is None):
                    resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1, _n_comp=numComp) #OC18072021
                    #resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
                if(resStokes3 is None):
                    resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1, _n_comp=numComp) #OC18072021
                    #resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
                #if(arAuxResSt123 is None):
                #    lenArSt123 = len(resStokes.arS) + len(resStokes2.arS) + len(resStokes3.arS)
                #    arAuxResSt123 = array('f', [0]*lenArSt123)
                lenArSt += (len(resStokes2.arS) + len(resStokes3.arS)) #OC24122018

            if(_pres_ang == 2): #24122018
                if((resStokesA is None) and (meshResA is not None)):
                    resStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny, _n_comp = numComp) #OC18072021
                    #resStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny)
                    lenArSt += len(resStokesA.arS) #OC26122018

            if((lenArSt > lenArSt0) and ((arAuxResSt is None) or (len(arAuxResSt) < lenArSt))): arAuxResSt = array('f', [0]*lenArSt) #OC24122018

            #DEBUG
            #print('workStokes.mesh: nx=', workStokes.mesh.nx, 'xStart=', workStokes.mesh.xStart, 'xFin=', workStokes.mesh.xFin)
            #print('workStokes.mesh: ny=', workStokes.mesh.ny, 'yStart=', workStokes.mesh.yStart, 'yFin=', workStokes.mesh.yFin)
            #print('resStokes.mesh: nx=', resStokes.mesh.nx, 'xStart=', resStokes.mesh.xStart, 'xFin=', resStokes.mesh.xFin)
            #print('resStokes.mesh: ny=', resStokes.mesh.ny, 'yStart=', resStokes.mesh.yStart, 'yFin=', resStokes.mesh.yFin)
            #END DEBUG
 
            if(_opt_bl is None):
                #resStokes.avg_update_same_mesh(workStokes, iAvgProc, 1)

                #print('resStokes.avg_update_same_mesh ... ', end='') #DEBUG
                #t0 = time.time(); #DEBUG
                
                if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021
                #if(_char != 6): #OC26022021
                    #resStokes.avg_update_same_mesh(workStokes, iAvgProc, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
                    resStokes.avg_update_same_mesh(workStokes, iAvgProc, numComp, ePhIntegMult) #OC16012017 #to treat all Stokes components / Polarization in the future

                    if((resStokes2 is not None) and (workStokes2 is not None)): #OC30052017
                        resStokes2.avg_update_same_mesh(workStokes2, iAvgProc, numComp, ePhIntegMult)

                    if((resStokes3 is not None) and (workStokes3 is not None)): #OC03052018
                        resStokes3.avg_update_same_mesh(workStokes3, iAvgProc, numComp, ePhIntegMult)

                    if((resStokesA is not None) and (workStokesA is not None)): #OC24122018
                        resStokesA.avg_update_same_mesh(workStokesA, iAvgProc, numComp, ePhIntegMult)

                #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                #DEBUG
                #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, 1)
                #END DEBUG
                
            else:
                #print('DEBUG MESSAGE: Started interpolation of current wavefront on resulting mesh')
                #if(doMutual <= 0): resStokes.avg_update_interp(workStokes, iAvgProc, 1, 1)
                #else: resStokes.avg_update_interp_mutual(workStokes, iAvgProc, 1)

                #print('resStokes.avg_update_interp ... ', end='') #DEBUG
                #t0 = time.time(); #DEBUG

                #DEBUG (test save at iAvgProc = 0)
                #if(iAvgProc == 0):
                #    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)
                #END DEBUG

                #if(doMutual <= 0): resStokes.avg_update_interp(workStokes, iAvgProc, 1, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future

                #if(_char == 40): #OC03052018
                if((_char == 40) or (_char == 41)): #OC15072019

                    #DEBUG
                    #if(i == 0):
                    #    print('Before first Stokes Averaging')
                    #    for ii in range(0, len(resStokes.arS), 10): print(ii, resStokes.arS[ii])

                    if(not (doPropCM and (iMode == iModeStart))): #OC23112022
                    #if(not (doPropCM and (iMode == 0))): #OC19112020
                    #if(not (doPropCM and (i == 0))): #OC13112020
                        
                        resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult, _sum=doPropCM) #OC04112020
                        #resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult)

                        if((resStokes2 is not None) and (workStokes2 is not None)):
                            #DEBUG
                            #print('Apdating from resStokes2')
                            resStokes2.avg_update_interp_mutual(workStokes2, iAvgProc, 1, ePhIntegMult, _sum=doPropCM) #OC13112020
                            #resStokes2.avg_update_interp_mutual(workStokes2, iAvgProc, 1, ePhIntegMult)

                        if((resStokes3 is not None) and (workStokes3 is not None)):
                            #DEBUG
                            #print('Apdating from resStokes3')
                            resStokes3.avg_update_interp_mutual(workStokes3, iAvgProc, 1, ePhIntegMult, _sum=doPropCM) #OC13112020
                            #resStokes3.avg_update_interp_mutual(workStokes3, iAvgProc, 1, ePhIntegMult)

                else:
                    if(doMutual == 0):
                        #DEBUG
                        #print('Before the update of Stokes')
                        #END DEBUG
                        #DEBUG
                        #print('Saving intensity of propagated wavefront #', i)
                        #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path + '_' + repr(i) + '.dat', 1, _mutual = doMutual)
                        #sys.exit()
                        #END DEBUG

                        if(not (doPropCM and (iMode == iModeStart))): #OC23112022
                        #if(not (doPropCM and (iMode == 0))): #OC19112020
                        #if(not (doPropCM and (i == 0))): #OC13112020

                            if resStokes.mesh.is_equal(workStokes.mesh): #OC01082022
                                resStokes.avg_update_same_mesh(workStokes, iAvgProc, numComp, ePhIntegMult, _sum=doPropCM)
                            else:
                                
                                try: #OC14012024 (added "try-except")
                                    resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult, _sum=doPropCM)
                                except:
                                    print('Failed to add intensity of single-electron / fully coherent radiation')

                            #resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult, _sum=doPropCM) #OC04112020
                        
                        elif resStokes.mesh.is_equal(workStokes.mesh): #OC01082022
                            resStokes.avg_update_same_mesh(workStokes, iAvgProc, numComp, ePhIntegMult, _sum=doPropCM)

                        #resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult) #OC16012017 #to treat all Stokes components / Polarization in the future
                        
                    else:
                        if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021 (added condition; MI was already updated in case of _char == 6)
                        #if(_char != 6): #OC04022021 (added condition; MI was already updated in case of _char == 6)
                            resStokes.avg_update_interp_mutual(workStokes, iAvgProc, 1, ePhIntegMult)

                        if((resStokes2 is not None) and (workStokes2 is not None)): #OC30052017
                            resStokes2.avg_update_interp_mutual(workStokes2, iAvgProc, 1, ePhIntegMult)

                if((resStokesA is not None) and (workStokesA is not None)): #OC24122018
                    resStokesA.avg_update_interp(workStokesA, iAvgProc, 1, numComp, ePhIntegMult)

                #DEBUG
                #srwl_uti_save_intens_ascii(resStokesA.arS, resStokesA.mesh, copy(_file_path) + '.ang_res.debug', 1)
                #END DEBUG
                
                #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                #print('DEBUG MESSAGE: Finished interpolation of current wavefront on resulting mesh')

            iAvgProc += 1
            #if(iAvgProc >= _n_part_avg_proc):

            doAllowSending = (nProc > 1) and (iAvgProc >= _n_part_avg_proc) #OC26022021
            #doAllowSending = (iAvgProc >= _n_part_avg_proc) #OC20112020
            if(doAllowSending and ((_char != 6) and (_char != 61) and (_char != 7))): #OC20062021 (consider removing the second part)
            #if(doAllowSending and (_char != 6)): #OC18022021 (consider removing the second part)
            #if(doAllowSending):
                if(doPropCM and ((nPartPerProc - i) <= _n_part_avg_proc) and ((nPartPerProc - i) > 1)): doAllowSending = False #OC20012021
                #if(doPropCM and ((nPartPerProc - i) < _n_part_avg_proc) and ((nPartPerProc - i) > 1)): doAllowSending = False #OC20112020

            #DEBUG
            #if(rank == 1): print('Rank #1: i=', i, 'iAvgProc=', iAvgProc, 'nPartPerProc=', nPartPerProc, '_n_part_avg_proc=', _n_part_avg_proc, 'doAllowSending=', doAllowSending)
            #sys.stdout.flush()
            #END DEBUG
            
            if(doAllowSending): #OC20112020
            #if(iAvgProc >= _n_part_avg_proc):

                #if(nProc > 1): #OC26022021: commented this out, since sending can be required only if nProc > 1

                #if(doPropCM and ((nPartPerProc - iMode) < _n_part_avg_proc) and ((nPartPerProc - iMode) > 1)): continue #OC20112020

                #sys.exit(0)
                #print("sending data from %d to 0" % rank) #an he
                #DEBUG
                #srwl_uti_save_intens_ascii(resStokes.arS, resStokes.mesh, _file_path, 1)
                #END DEBUG

                #DEBUG
                #srwl_uti_save_text("Preparing to sending # " + str(iAuxSendCount + 1), _file_path + "." + str(rank) + "bs.dbg")
                #END DEBUG

                ##comMPI.Send([resStokes.arS, MPI.FLOAT], dest=0)
                #if(resStokes2 is None):
                #    comMPI.Send([resStokes.arS, MPI.FLOAT], dest=0)
                ##else: #OC31052017
                #elif(resStokes3 is None): #OC03052018
                #    lenArSt1 = len(resStokes.arS)
                #    lenArSt2 = len(resStokes2.arS)
                #    for i1 in range(lenArSt1): arAuxResSt12[i1] = resStokes.arS[i1]
                #    for i2 in range(lenArSt2): arAuxResSt12[i2 + lenArSt1] = resStokes2.arS[i2]
                #    comMPI.Send([arAuxResSt12, MPI.FLOAT], dest=0)
                #else: #OC03052018
                #    lenArSt1 = len(resStokes.arS)
                #    lenArSt2 = len(resStokes2.arS)
                #    lenArSt3 = len(resStokes3.arS)
                #    for i1 in range(lenArSt1): arAuxResSt123[i1] = resStokes.arS[i1]
                #    for i2 in range(lenArSt2): arAuxResSt123[i2 + lenArSt1] = resStokes2.arS[i2]
                #    lenArSt12 = lenArSt1 + lenArSt2
                #    for i3 in range(lenArSt3): arAuxResSt123[i3 + lenArSt12] = resStokes3.arS[i3]
                #    comMPI.Send([arAuxResSt123, MPI.FLOAT], dest=0)

                if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021
                #if(_char != 6): #OC18022021 (in case _char == 6 Electric Field should be sent)
                    #OC24122018
                    resStkToSend = resStokes.arS
                    lenArSt1 = len(resStokes.arS)
                    lenArSt = lenArSt1
                    if((resStokes2 is not None) and (resStokes3 is None)):
                        lenArSt2 = len(resStokes2.arS)
                        #for i1 in range(lenArSt): arAuxResSt[i1] = resStokes.arS[i1]
                        arAuxResSt[0:lenArSt] = resStokes.arS #??
                        #for i2 in range(lenArSt2): arAuxResSt[i2 + lenArSt] = resStokes2.arS[i2]
                        lenArStTot = lenArSt + lenArSt2
                        arAuxResSt[lenArSt:lenArStTot] = resStokes2.arS #??
                        lenArSt = lenArStTot
                        resStkToSend = arAuxResSt
                    elif((resStokes2 is not None) and (resStokes3 is not None)):
                        lenArSt2 = len(resStokes2.arS)
                        lenArSt3 = len(resStokes3.arS)
                        #for i1 in range(lenArSt): arAuxResSt[i1] = resStokes.arS[i1]
                        arAuxResSt[0:lenArSt] = resStokes.arS #??
                        #for i2 in range(lenArSt2): arAuxResSt[i2 + lenArSt] = resStokes2.arS[i2]
                        lenArStTot = lenArSt + lenArSt2
                        arAuxResSt[lenArSt:lenArStTot] = resStokes2.arS #??
                        lenArSt = lenArStTot
                        #for i3 in range(lenArSt3): arAuxResSt[i3 + lenArSt] = resStokes3.arS[i3]
                        lenArStTot = lenArSt + lenArSt3
                        arAuxResSt[lenArSt:lenArStTot] = resStokes3.arS #??
                        lenArSt = lenArStTot
                        resStkToSend = arAuxResSt

                    if(resStokesA is not None):
                        if(resStokes2 is None):
                            #for i1 in range(lenArSt): arAuxResSt[i1] = resStokes.arS[i1]
                            arAuxResSt[0:lenArSt] = resStokes.arS
                        lenArStA = len(resStokesA.arS)
                        #for ia in range(lenArStA): arAuxResSt[ia + lenArSt] = resStokesA.arS[ia]
                        lenArStTot = lenArSt + lenArStA
                        arAuxResSt[lenArSt:lenArStTot] = resStokesA.arS
                        lenArSt = lenArStTot
                        resStkToSend = arAuxResSt
                        
                    #comMPI.Send([arAuxResSt, MPI.FLOAT], dest=0)
                    #comMPI.Send([resStkToSend, MPI.FLOAT], dest=0) #OC26122018
                    comMPI.Send([resStkToSend, MPI.FLOAT], dest=rankMaster) #OC03032021

                    #if(resStokes2 is not None): comMPI.Send([resStokes2.arS, MPI.FLOAT], dest=0) #OC30052017

                    for ir in range(len(resStokes.arS)): #OC17022021
                    #for ir in range(nStPt):
                        resStokes.arS[ir] = 0

                    if(resStokes2 is not None): #OC30052017
                        for ir in range(len(resStokes2.arS)): #OC17022021
                        #for ir in range(nStPt2):
                            resStokes2.arS[ir] = 0

                    if(resStokes3 is not None): #OC03052018
                        for ir in range(len(resStokes3.arS)): #OC17022021
                        #for ir in range(nStPt3):
                            resStokes3.arS[ir] = 0

                    if(resStokesA is not None): #OC27122018
                        for ir in range(len(resStokesA.arS)): resStokesA.arS[ir] = 0

                else: #if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021 (in case _char == 6 Electric Field should be sent)
                #else: #if(_char == 6): #OC18022021 (in case _char == 6 Electric Field should be sent)

                    #Resize the Electric Field according to the required meshRes (commented-out, since this was done before)
                    #srwl.ResizeElecFieldMesh(wfr, meshRes, [0,1]) #[0,1] means do the resizing without FFT and allow treatment of quad. phase terms

                    #DEBUG
                    #print('Rank #', rank, ': mode #', iMode, 'about to send electric field calculated on the mesh:')
                    #print('wfr.mesh.eStart=', wfr.mesh.eStart, 'wfr.mesh.eFin=', wfr.mesh.eFin, 'wfr.mesh.ne=', wfr.mesh.ne)
                    #print('wfr.mesh.xStart=', wfr.mesh.xStart, 'wfr.mesh.xFin=', wfr.mesh.xFin, 'wfr.mesh.nx=', wfr.mesh.nx)
                    #print('wfr.mesh.yStart=', wfr.mesh.yStart, 'wfr.mesh.yFin=', wfr.mesh.yFin, 'wfr.mesh.ny=', wfr.mesh.ny)
                    #sys.stdout.flush()
                    #t0 = time.time()
                    #END DEBUG

                    #Consider adding logic if wfr.arEx or wfr.arEy is not defined
                    arElFldToSend[0:lenHalfArToSend] = wfr.arEx
                    arElFldToSend[lenHalfArToSend:lenArToSend] = wfr.arEy

                    #DEBUG_OC16042021 (commented-out comMPI.Send([arElFldToSend, MPI.FLOAT], dest=rankMaster) line for test)
                    comMPI.Send([arElFldToSend, MPI.FLOAT], dest=rankMaster) #OC03032021 (sending electric field to master)
                    #comMPI.Send([arElFldToSend, MPI.FLOAT], dest=0) #OC18022021 (sending electric field to master)

                    #DEBUG
                    #print('Rank #', rank, ': mode #', iMode, 'sent to Master; sending lasted:', round(time.time() - t0, 6), 's')
                    #sys.stdout.flush()
                    #END DEBUG

                #OC18022021 (moved here from upper location)
                iAuxSendCount += 1 #for debug

                #DEBUG
                #srwl_uti_save_text("Sent # " + str(iAuxSendCount), _file_path + "." + str(rank) + "es.dbg")
                #END DEBUG

                #DEBUG
                #srwl_uti_save_intens_ascii(resStokes.arS, resStokes.mesh, _file_path, 1)
                #print('Rank #', rank, ': summed-up propagated mode data sent to Master; i=', i, ' iAvgProc=', iAvgProc)
                #sys.stdout.flush()
                #END DEBUG

                iAvgProc = 0

            if(nProc == 1):
                #DEBUG
                #if(i == 1):
                #    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1)
                #    sys.exit(0)
                #END DEBUG
                iSave += 1
                if((_file_path is not None) and (iSave == _n_save_per)):
                    #Saving results from time to time in the process of calculation:

                    #print('srwl_uti_save_intens_ascii ... ', end='') #DEBUG
                    #t0 = time.time(); #DEBUG
                    
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC26042016
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017
                        
                    fp = _file_path; fp1 = file_path1; fp2 = file_path2; fpdc1 = file_path_deg_coh1; fpdc2 = file_path_deg_coh2 #OC14082018
                    fpA = file_pathA #OC24122018
                    if(_file_bkp): 
                        if(bkpFileToBeSaved):
                            if(fp is not None): fp = copy(fp) + '.bkp'
                            if(fp1 is not None): fp1 = copy(fp1) + '.bkp'
                            if(fp2 is not None): fp2 = copy(fp2) + '.bkp'
                            if(fpdc1 is not None): fpdc1 = copy(fpdc1) + '.bkp'
                            if(fpdc2 is not None): fpdc2 = copy(fpdc2) + '.bkp'
                            if(fpA is not None): fpA = copy(fpA) + '.bkp' #OC24122018
                            
                            bkpFileToBeSaved = False
                        else: bkpFileToBeSaved = True

                    if(((_char == 6) or (_char == 61) or (_char == 7)) and (_n_mpi <= 1)): #OC20062021 (copy / update CSD only if total distribution is required)
                    #if((_char == 6) and (_n_mpi <= 1)): #OC03032021 (copy / update CSD only if total distribution is required)
                    #if(_char == 6):
                        #DEBUG
                        #print('About to fill symmetrical part of Hermitian Mutual Intensity distribution')
                        #END DEBUG

                        srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [4]) #Filling-in "symmetrical" part of the Hermitian Mutual Intensity distribution

                    #if(_char == 40): #OC03052018
                    if((_char == 40) or (_char == 41)): #OC13072019
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0)
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, fp, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0) #OC14082018 #Intensity
                        srwl_uti_save_intens(resStokes.arS, meshRes, fp, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0, _form = _file_form) #OC17072021 #Intensity

                        if(resStokes2 is not None):
                            #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                            #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                            if(_char == 40): #OC13072019
                                srwl_uti_save_intens(resStokes2.arS, resStokes2.mesh, fp1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1, _form = _file_form) #OC17072021 #Mutual Intensity, Hor. Cut
                                #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, fp1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC14082018 #Mutual Intensity, Hor. Cut
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC06052018
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC14082018 #Degree of Coherence, Hor. Cut
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 #Degree of Coherence, Hor. Cut

                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 #Degree of Coherence, Hor. Cut
                            #OCTEST14112020
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(_rel_zer_tol=0), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 #Degree of Coherence, Hor. Cut
                            srwl_uti_save_intens(resStokes2.to_deg_coh(_rel_zer_tol=0), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 #Degree of Coherence, Hor. Cut

                        if(resStokes3 is not None):
                            #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                            #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                            if(_char == 40): #OC13072019
                                srwl_uti_save_intens(resStokes3.arS, resStokes3.mesh, fp2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1, _form = _file_form) #OC17072021 #Mutual Intensity, Vert. Cut
                                #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, fp2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC14082018 #Mutual Intensity, Vert. Cut
 
                            #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC06052018
                            #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC14082018 #Degree of Coherence, Vert. Cut
                            #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 #Degree of Coherence, Vert. Cut

                            #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 #Degree of Coherence, Vert. Cut
                            #OCTEST14112020
                            #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(_rel_zer_tol=0), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 #Degree of Coherence, Vert. Cut
                            srwl_uti_save_intens(resStokes3.to_deg_coh(_rel_zer_tol=0), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 #Degree of Coherence, Vert. Cut

                    #elif(_char == 4): #OC03052018
                    elif((_char == 4) or (_char == 5)): #OC13072019
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                        if(_char == 4):
                            srwl_uti_save_intens(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021 #Mutual Intensity, Hor. Cut
                            #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018 #Mutual Intensity, Hor. Cut
                        #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC06052018
                        #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC14082018 #Degree of Coherence, Hor. Cut
                        #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 #Degree of Coherence, Hor. Cut
                        #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 #Degree of Coherence, Hor. Cut
                        srwl_uti_save_intens(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 #Degree of Coherence, Hor. Cut
                        if((resStokes2 is not None) and (meshRes2 is not None)): #OC30052017
                            #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                            #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                            if(_char == 4):
                                srwl_uti_save_intens(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021  #Mutual Intensity, Vert. Cut
                                #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018  #Mutual Intensity, Vert. Cut
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0) #OC06052018
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0) #OC14082018 #Degree of Coherence, Vert. Cut
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = doMutual, _cmplx = 0) #OC12072019 #Degree of Coherence, Vert. Cut
                            #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 #Degree of Coherence, Vert. Cut
                            srwl_uti_save_intens(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 #Degree of Coherence, Vert. Cut

                    else:
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                        #DEBUG
                        #print('meshRes.eStart=', meshRes.eStart, ' meshRes.eFin=', meshRes.eFin)
                        #END DEBUG

                        srwl_uti_save_intens(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form, _wfr = wfrA) #OC18062021
                        #srwl_uti_save_intens(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC18022021
                        #if(_file_form == 'ascii'): #OC05022021
                        #    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                        #elif(_file_form == 'hdf5'): #OC05022021
                        #    srwl_uti_save_intens_hdf5(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                        
                        if((resStokes2 is not None) and (meshRes2 is not None)): #OC30052017
                            #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                            #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                            srwl_uti_save_intens(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021

                    if(_pres_ang == 2): #OC24122018
                        srwl_uti_save_intens(resStokesA.arS, meshResA, fpA, numComp, _arLabels = resLabelsToSaveA, _arUnits = resUnitsToSaveA, _mutual = 0, _cmplx = 0, _form = _file_form) #OC17072021
                        #srwl_uti_save_intens_ascii(resStokesA.arS, meshResA, fpA, numComp, _arLabels = resLabelsToSaveA, _arUnits = resUnitsToSaveA, _mutual = 0, _cmplx = 0)
                        
                    #DEBUG
                    #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)
                    #END DEBUG

                    #MR01112016: write the status of the simulation:  
                    #srwl_uti_save_stat_wfr_emit_prop_multi_e(i + 1, total_num_of_particles, filename=log_path)
                    srwl_uti_save_stat_wfr_emit_prop_multi_e(i + 1, actNumPartTot, filename=log_path) #OC31102021
                    
                    #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                    
                    #sys.exit(0)
                    iSave = 0

    elif((rank == rankMaster) and (nProc > 1)): #OC02032021
    #elif((rank == 0) and (nProc > 1)):

        #nRecv = int(nPartPerProc*nProc/_n_part_avg_proc + 1e-09)
        nRecv = nSentPerProc*(nProc - 1) #Total number of sending acts to be made by all worker processes, and to be received by master

        if(doPropCM and ((_char == 6) or (_char == 61) or (_char == 7))): nRecv = _n_part_tot #OC28102021

        #print('DEBUG MESSAGE: Actual number of macro-electrons:', nRecv*_n_part_avg_proc)
        
        #DEBUG
        #srwl_uti_save_text("nRecv: " + str(nRecv) + " nPartPerProc: " + str(nPartPerProc) + " nProc: " + str(nProc) + " _n_part_avg_proc: " + str(_n_part_avg_proc), _file_path + ".00.dbg")

        #print('Master: nRecv=', nRecv)
        #sys.stdout.flush()
        #END DEBUG

        if(resStokes is None):
            #nComp = 4 #OC19022021
            #if((_char == 6) or (_char == 61) or (_char == 7)): nComp = 1 #OC20062021
            ##if(_char == 6): nComp = 1 #OC19022021
            #OC18072021 (commenyed-out the above, using numComp instead)

            #OC06052023
            goAlloc = False
            if((_char == 6) or (_char == 61) or (_char == 7)):
                if(rank == 0): goAlloc = True
                elif(rank == rankMaster):
                    auxRecv = array('i', [0])
                    comMPI.Recv([auxRecv, MPI.INT], source=(rank - nProc))
                    if(auxRecv[0] > 0): goAlloc = True
            else:
                goAlloc = True
            
            if(goAlloc): #OC06052023
                resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, _n_comp = numComp, _itStFin=itStartEnd)
                if(rank < (nProcTot - nProc)):
                    auxSend = array('i', [1])
                    comMPI.Send([auxSend, MPI.INT], dest=(rank + nProc))
                    
                #DEBUG
                #print('rank=', rank, 'memory for CSD allocated')
                #sys.stdout.flush()
                #END DEBUG
            
            #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, _n_comp = nComp, _itStFin=itStartEnd) #OC03032021
            #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, _n_comp = nComp) #OC19022021
            #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

            #DEBUG
            #print('Rank #', rank, ': resStokes allocated:')
            #print('meshRes.eStart=', meshRes.eStart, 'meshRes.eFin=', meshRes.eFin, 'meshRes.ne=', meshRes.ne)
            #print('meshRes.xStart=', meshRes.xStart, 'meshRes.xFin=', meshRes.xFin, 'meshRes.nx=', meshRes.nx)
            #print('meshRes.yStart=', meshRes.yStart, 'meshRes.yFin=', meshRes.yFin, 'meshRes.ny=', meshRes.ny)
            #sys.stdout.flush()
            #END DEBUG

        if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021
        #if(_char != 6): #OC18022021
            if(workStokes is None):
                workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, _n_comp = numComp) #OC18072021
                #workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

        if((_char == 6) or (_char == 61) or (_char == 7)): #OC03102021
        #if(((_char == 6) or (_char == 61) or (_char == 7)) and (_opt_bl is None)): #OC20062021
        #if((_char == 6) and (_opt_bl is None)): #OC18062021

            if(wfrA is None): wfrA = SRWLWfr()

            if(wfrA.avgPhotEn <= 0): #I.e. if Average Wavefront was not filled-out

                if doPropCM: #OC26102021
                    wfr = _mag[0]

                    #OC10092022
                    if(_e_beam is not None):
                        m1 = _e_beam.partStatMom1; m1w = wfr.partBeam.partStatMom1
                        dx = m1.x - m1w.x; dxp = m1.xp - m1w.xp; dy = m1.y - m1w.y; dyp = m1.yp - m1w.yp
                        if((dx != 0.) or (dxp != 0.) or (dy != 0.) or (dyp != 0.)): 
                            wfr.sim_src_offset(_dx = dx, _dxp = dxp, _dy = dy, _dyp = dyp, _move_mesh=False, _copy=False)
                    #OC01082022
                    #if((wfr.partBeam.partStatMom1.x != _e_beam.partStatMom1.x) or (wfr.partBeam.partStatMom1.xp != _e_beam.partStatMom1.xp) or 
                    #   (wfr.partBeam.partStatMom1.y != _e_beam.partStatMom1.y) or (wfr.partBeam.partStatMom1.yp != _e_beam.partStatMom1.yp)):
                    #    wfr.sim_src_offset(_dx = (_e_beam.partStatMom1.x - wfr.partBeam.partStatMom1.x), _dxp = (_e_beam.partStatMom1.xp - wfr.partBeam.partStatMom1.xp), 
                    #                       _dy = (_e_beam.partStatMom1.y - wfr.partBeam.partStatMom1.y), _dyp = (_e_beam.partStatMom1.yp - wfr.partBeam.partStatMom1.yp), _move_mesh=False, _copy=False)
                        #wfr = wfr.sim_src_offset(_dx = (_e_beam.partStatMom1.x - wfr.partBeam.partStatMom1.x), _dxp = (_e_beam.partStatMom1.xp - wfr.partBeam.partStatMom1.xp), 
                        #                         _dy = (_e_beam.partStatMom1.y - wfr.partBeam.partStatMom1.y), _dyp = (_e_beam.partStatMom1.yp - wfr.partBeam.partStatMom1.yp), _move_mesh=False, _copy=True)
                else:
                    srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR) #It's unfortunate to calculate this here, but we need to know Rx, Ry, etc.

                if(_opt_bl is not None): srwl.PropagElecField(wfr, _opt_bl) #OC05102021 (for cases when CMD has to be done after propagation)
                #OC05102021: furthemore, wfr.arEx, wfr.arEy are necessary for receiving E-field data from workers and accumulating CSD from it

                if(_det is not None): srwl.ResizeElecFieldMesh(wfr, meshRes, [0, 1]) #OC05102021

                wfrA.Rx = wfr.Rx
                wfrA.Ry = wfr.Ry
                wfrA.dRx = wfr.dRx
                wfrA.dRy = wfr.dRy

                if doPropCM: #OC26102021
                    wfrA.xc = wfr.xc
                    wfrA.yc = wfr.yc
                else:
                    wfrA.xc = elecX0
                    wfrA.yc = elecY0

                wfrA.avgPhotEn = wfr.avgPhotEn
                wfrA.presCA = wfr.presCA
                wfrA.presFT = wfr.presFT
                wfrA.unitElFld = wfr.unitElFld
                wfrA.unitElFldAng = wfr.unitElFldAng
                wfrA.mesh = copy(wfr.mesh) #OC27102021
                #wfrA.mesh = wfr.mesh #OC28062021

        lenArResSt0 = len(resStokes.arS) #OC24122018
        lenArResSt = lenArResSt0
        
        lenArWorkSt0 = 0 #OC18022021
        if(workStokes is not None): #OC18022021
            lenArWorkSt0 = len(workStokes.arS)
        lenArWorkSt = lenArWorkSt0

        #OC18022021
        arElFldToRecv = None #To be used only in the case on nProc > 1 and _char == 6
        lenHalfArToRecv = 0; lenArToRecv = 0
        if((_char == 6) or (_char == 61) or (_char == 7)): #OC20062021
        #if(_char == 6): #To allocate array for sending Electric Field, using the numbers of points in the final mesh (to be taken from _det or from first test propagated wavefront)
            #Asumes that meshRes is already defined at this point for workers
            lenHalfArToRecv = meshRes.ne*meshRes.nx*meshRes.ny*2 #for arEx and arEy (consider introducing some logic of arEx or arEy is not used)
            lenArToRecv = 2*lenHalfArToRecv
            arElFldToRecv = array('f', [0]*lenArToRecv)
        
        #if(_char == 4): #OC30052017 #Cuts of Mutual Intensity vs X & Y
        if((_char == 4) or (_char == 5)): #OC15072019 #Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
            if(resStokes2 is None):
                resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual, _n_comp = numComp) #OC18072021
                #resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
            #if(arAuxResSt12 is None):
            #    lenArResSt12 = len(resStokes.arS) + len(resStokes2.arS)
            #    arAuxResSt12 = array('f', [0]*lenArResSt12)
            lenArResSt += len(resStokes2.arS) #OC24122028

            if(workStokes2 is None):
                workStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual, _n_comp = numComp) #OC18072021
                #workStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
            #if(arAuxWorkSt12 is None):
            #    lenArWorkSt12 = len(workStokes.arS) + len(workStokes2.arS)
            #    arAuxWorkSt12 = array('f', [0]*lenArWorkSt12)
            lenArWorkSt += len(workStokes2.arS) #OC24122018

        #if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
        if((_char == 40) or (_char == 41)): #OC15072019 #Intensity and Cuts of Mutual Intensity and/or Degree of Coherence vs X & Y
            if(resStokes2 is None):
                resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1, _n_comp = numComp) #OC18072021
                #resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
            if(resStokes3 is None):
                resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1, _n_comp = numComp) #OC18072021
                #resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
            #if(arAuxResSt123 is None):
            #    lenArResSt123 = len(resStokes.arS) + len(resStokes2.arS) + len(resStokes3.arS)
            #    arAuxResSt123 = array('f', [0]*lenArResSt123)
            lenArResSt += (len(resStokes2.arS) + len(resStokes3.arS)) #OC24122018

            if(workStokes2 is None):
                workStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1, _n_comp = numComp) #OC18072021
                #workStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
            if(workStokes3 is None):
                workStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1, _n_comp = numComp) #OC18072021
                #workStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
            #if(arAuxWorkSt123 is None):
            #    lenArWorkSt123 = len(workStokes.arS) + len(workStokes2.arS) + len(workStokes3.arS)
            #    arAuxWorkSt123 = array('f', [0]*lenArWorkSt123)
            lenArWorkSt += (len(workStokes2.arS) + len(workStokes3.arS)) #OC24122018

        if(_pres_ang == 2): #OC24122018
            if(resStokesA is None):
                resStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny, _n_comp = numComp) #OC18072021
                #resStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny)
            if(workStokesA is None):
                workStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny, _n_comp = numComp) #OC18072021
                #workStokesA = SRWLStokes(1, 'f', meshResA.eStart, meshResA.eFin, meshResA.ne, meshResA.xStart, meshResA.xFin, meshResA.nx, meshResA.yStart, meshResA.yFin, meshResA.ny)
            lenArResSt += len(resStokesA.arS)
            lenArWorkSt += len(workStokesA.arS)

        workStkToRcv = None #OC18022021
        if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021
        #if(_char != 6): #OC18022021
            workStkToRcv = workStokes.arS #OC24122018

            if(lenArResSt > lenArResSt0): arAuxResSt = array('f', [0]*lenArResSt) #OC24122018
            if(lenArWorkSt > lenArWorkSt0):
                arAuxWorkSt = array('f', [0]*lenArWorkSt)
                workStkToRcv = arAuxWorkSt

        else: #if((_char == 6) or (_char == 61) or (_char == 7)) #OC20062021 (allocating aux. wrf to be able to extract CSD from received Electric Field data)
        #else: #if(_char == 6) #OC19022021 (allocating aux. wrf to be able to extract CSD from received Electric Field data)
            if(wfr is None): #OC16042021 (maybe the thing below is not at all necessary?)
                wfr = SRWLWfr(_arEx=1, _arEy=1, _typeE='f', _eStart=meshRes.eStart, _eFin=meshRes.eFin, _ne=meshRes.ne, 
                              _xStart=meshRes.xStart, _xFin=meshRes.xFin, _nx=meshRes.nx, 
                              _yStart=meshRes.yStart, _yFin=meshRes.yFin, _ny=meshRes.ny)

        for i in range(nRecv): #loop over messages from workers

            #DEBUG
            #srwl_uti_save_text("Preparing to receiving # " + str(i), _file_path + ".br.dbg")
            #END DEBUG
           
            ##comMPI.Recv([workStokes.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 
            #if(workStokes2 is None): #OC31052017
            #    comMPI.Recv([workStokes.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 
            ##else:
            #elif(workStokes3 is None): #OC03052018
            #    comMPI.Recv([arAuxWorkSt12, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive
            #
            #    lenArWorkSt1 = len(workStokes.arS)
            #    for i1 in range(lenArWorkSt1): workStokes.arS[i1] = arAuxWorkSt12[i1]
            #    lenArWorkSt2 = len(workStokes2.arS)
            #    for i2 in range(lenArWorkSt2): workStokes2.arS[i2] = arAuxWorkSt12[i2 + lenArWorkSt1]
            #
            #else: #OC03052018
            #    comMPI.Recv([arAuxWorkSt123, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 
            #
            #    lenArWorkSt1 = len(workStokes.arS)
            #    for i1 in range(lenArWorkSt1): workStokes.arS[i1] = arAuxWorkSt123[i1]
            #    lenArWorkSt2 = len(workStokes2.arS)
            #    for i2 in range(lenArWorkSt2): workStokes2.arS[i2] = arAuxWorkSt123[i2 + lenArWorkSt1]
            #    lenArWorkSt3 = len(workStokes3.arS)
            #    lenArWorkSt12 = lenArWorkSt1 + lenArWorkSt2
            #    for i3 in range(lenArWorkSt3): workStokes3.arS[i3] = arAuxWorkSt123[i3 + lenArWorkSt12]

            if((_char != 6) and (_char != 61) and (_char != 7)): #OC20062021
            #if(_char != 6): #OC18022021

                #DEBUG
                #print('About to start receiving by Master at i=', i)
                #sys.stdout.flush()
                #END DEBUG
                #OC24122018
                comMPI.Recv([workStkToRcv, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 

                lenArSt1 = len(workStokes.arS)
                lenArSt = lenArSt1

                if((workStokes2 is not None) and (workStokes3 is None)):
                    for i1 in range(lenArSt): workStokes.arS[i1] = arAuxWorkSt[i1]
                    #workStokes.arS[0:lenArSt] = arAuxWorkSt
                    lenArSt2 = len(workStokes2.arS)
                    for i2 in range(lenArSt2): workStokes2.arS[i2] = arAuxWorkSt[i2 + lenArSt]
                    lenArSt += lenArSt2

                elif((workStokes2 is not None) and (workStokes3 is not None)):
                    for i1 in range(lenArSt): workStokes.arS[i1] = arAuxWorkSt[i1]
                    #workStokes.arS[0:lenArSt] = arAuxWorkSt
                    lenArSt2 = len(workStokes2.arS)
                    for i2 in range(lenArSt2): workStokes2.arS[i2] = arAuxWorkSt[i2 + lenArSt]
                    lenArSt += lenArSt2
                    lenArSt3 = len(workStokes3.arS)
                    for i3 in range(lenArSt3): workStokes3.arS[i3] = arAuxWorkSt[i3 + lenArSt]
                    lenArSt += lenArSt3

                if(workStokesA is not None): #OC24122018
                    #workStokes.arS[0:lenArSt] = arAuxWorkSt
                    if(workStokes2 is None): #OC26122018
                        for j in range(lenArSt): workStokes.arS[j] = arAuxWorkSt[j] #OC26122018
                
                    lenArStA = len(workStokesA.arS)
                    for ia in range(lenArStA): workStokesA.arS[ia] = arAuxWorkSt[ia + lenArSt]
                    lenArSt += lenArStA

                #DEBUG
                #print('Propagated mode data received by Master. Number of modes received so far:', i+1)
                #sys.stdout.flush()
                #END DEBUG

                #if((_char == 4) and (workStokes2 is not None)): #OC30052017 #Cuts of Mutual Intensity vs X & Y
                #    comMPI.Recv([workStokes2.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 

                #OC15102018 (moved this log-writing to the place where other files are saved)
                #MR20160907 #Save .log and .json files:
                #particle_number = (i + 1) * _n_part_avg_proc
                #srwl_uti_save_stat_wfr_emit_prop_multi_e(particle_number, total_num_of_particles, filename=log_path)

                #DEBUG
                #srwl_uti_save_text("Received intensity # " + str(i), _file_path + ".er.dbg")
                #END DEBUG

                #resStokes.avg_update_same_mesh(workStokes, i + 1)
                #resStokes.avg_update_same_mesh(workStokes, i + 1, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
                multFinAvg = 1 if(_n_part_avg_proc > 1) else ePhIntegMult #OC120714 fixed: the normalization may have been already applied at the previous averaging in each worker process!

                #print('resStokes.avg_update_same_mesh ... ', end='') #DEBUG
                #t0 = time.time(); #DEBUG

                #DEBUG (test save at i = 0)
                #if(i == 0):
                #    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path + '.0.dat', 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)
                #END DEBUG

                #resStokes.avg_update_same_mesh(workStokes, i + 1, 1, multFinAvg) #in the future treat all Stokes components / Polarization, not just s0!
                #resStokes.avg_update_same_mesh(workStokes, i, 1, multFinAvg) #OC15012017 #in the future treat all Stokes components / Polarization, not just s0!
                #resStokes.avg_update_same_mesh(workStokes, i, numComp, multFinAvg) #OC15012017 #in the future treat all Stokes components / Polarization, not just s0!
                resStokes.avg_update_same_mesh(workStokes, i, numComp, multFinAvg, _sum=doPropCM) #OC20112020

                if((resStokes2 is not None) and (workStokes2 is not None)):
                    resStokes2.avg_update_same_mesh(workStokes2, i, numComp, multFinAvg, _sum=doPropCM) #OC20112020
                    #resStokes2.avg_update_same_mesh(workStokes2, i, numComp, multFinAvg) #OC30052017 #in the future treat all Stokes components / Polarization, not just s0!

                if((resStokes3 is not None) and (workStokes3 is not None)):
                    resStokes3.avg_update_same_mesh(workStokes3, i, numComp, multFinAvg, _sum=doPropCM) #OC20112020
                    #resStokes3.avg_update_same_mesh(workStokes3, i, numComp, multFinAvg) #OC03052018 #in the future treat all Stokes components / Polarization, not just s0!

                if((resStokesA is not None) and (workStokesA is not None)): #OC24122018
                    resStokesA.avg_update_same_mesh(workStokesA, i, numComp, multFinAvg, _sum=doPropCM) #OC20112020
                    #resStokesA.avg_update_same_mesh(workStokesA, i, numComp, multFinAvg) #in the future treat all Stokes components / Polarization, not just s0!

                #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                #DEBUG
                #srwl_uti_save_text("Updated Stokes after receiving intensity # " + str(i), _file_path + "." + str(i) + "er.dbg")
                #END DEBUG

            else: #OC20062021 if((_char == 6) or (_char == 61) or (_char == 7))
            #else: #OC18022021 if(_char == 6)

                #DEBUG
                #import numpy as np
                #print('rank=', rank, 'iRecv=', i, ': Trying to receive E-field')
                #sys.stdout.flush()
                #END DEBUG

                #DEBUG_OC16042021: commented-out the line: comMPI.Recv([arElFldToRecv, MPI.FLOAT], source=MPI.ANY_SOURCE) for testing
                comMPI.Recv([arElFldToRecv, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive Electric Field

                #DEBUG
                #print('rank=', rank, ': Electric Field data received by Master. Number of Electric Fields received so far:', i+1)
                #sys.stdout.flush()
                #END DEBUG

                #OC18022021
                wfr.arEx[:] = arElFldToRecv[0:lenHalfArToRecv]
                wfr.arEy[:] = arElFldToRecv[lenHalfArToRecv:lenArToRecv]
                #Any other params need to be set in wfr?

                #DEBUG OC16102021
                #import numpy as np
                #print('Trying to test recieve E-field iRecv=', i)
                #np.asarray_chkfinite(wfr.arEx, dtype=float)
                #np.asarray_chkfinite(wfr.arEy, dtype=float)
                #print('rank=', rank, 'iRecv=', i, ': No NaN or Inf received in E-field')
                #sys.stdout.flush()
                #END DEBUG

                #Update 4D CSD from the received Electric Field data
                intSumType = 1 #calculation of Intensity with instant averaging
                if(doPropCM): intSumType = 2 #adding of new Intensity value to previous one

                #DEBUG
                #t0 = time.time()
                #END DEBUG

                arMethPar = [0]*20 #OC03032021 (in principle, this is not required if(nProc == 1) ?)
                arMethPar[0] = intSumType; arMethPar[1] = i
                if((_n_mpi > 1) and (itStartEnd is not None)):
                    arMethPar[18] = itStartEnd[0]
                    arMethPar[19] = itStartEnd[1]

                srwl.CalcIntFromElecField(resStokes.arS, wfr, -1, 8, depTypeInt, phEnInt, 0., 0., arMethPar) #OC03032021
                #srwl.CalcIntFromElecField(resStokes.arS, wfr, 6, 8, depTypeInt, phEnInt, 0., 0., [intSumType, i]) #One main Stokes component

                #DEBUG OC16102021
                #import numpy as np
                #try:
                #    np.asarray_chkfinite(resStokes.arS, dtype=float)
                #except:
                #    traceback.print_exc()
                #    for iii in range(len(resStokes.arS)):
                #        if(isnan(resStokes.arS[iii])): 
                #            print('NaN found in resStokes.arS at iii=', iii)
                #            break
                #        elif(isinf(resStokes.arS[iii])): 
                #            print('Inf found in resStokes.arS at iii=', iii)
                #            break
                #print('rank=', rank, 'iRecv=', i, ': No NaN or Inf found in resStokes.arS (instant CSD)')
                #sys.stdout.flush()
                #END DEBUG

                #DEBUG
                #print('rank=', rank, ': Propagated mode data received by Master and CSD updated. Number of Electric Fields received so far:', i, 'Update lasted:', round(time.time() - t0, 6), 's')
                #sys.stdout.flush()
                #END DEBUG

            iSave += 1
            if(iSave == _n_save_per):
                #Saving results from time to time in the process of calculation

                #print('srwl_uti_save_intens_ascii ... ', end='') #DEBUG
                #t0 = time.time(); #DEBUG
                
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC26042016
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017

                #OC15102018 (moved this log-writing to the place where other files are saved):
                #MR20160907 #Save .log and .json files:
                particle_number = (i + 1) * _n_part_avg_proc
                if i == (nRecv - 1): particle_number = actNumPartTot #OC31102021
                #if i == (nRecv - 1): particle_number = total_num_of_particles #OC21012021 (previous particle_numbers may be inaccurate)
                srwl_uti_save_stat_wfr_emit_prop_multi_e(particle_number, actNumPartTot, filename=log_path) #OC31102021
                #srwl_uti_save_stat_wfr_emit_prop_multi_e(particle_number, total_num_of_particles, filename=log_path)

                #DEBUG
                #print('Updated log-file: particle_number=', particle_number, ' actNumPartTot=', actNumPartTot)
                #sys.stdout.flush()
                #END DEBUG

                fp = _file_path; fp1 = file_path1; fp2 = file_path2; fpdc1 = file_path_deg_coh1; fpdc2 = file_path_deg_coh2 #OC14082018
                fpA = file_pathA #OC24122018
                if(_file_bkp): 
                    if(bkpFileToBeSaved):
                        if(fp is not None): fp = copy(fp) + '.bkp'
                        if(fp1 is not None): fp1 = copy(fp1) + '.bkp'
                        if(fp2 is not None): fp2 = copy(fp2) + '.bkp'
                        if(fpdc1 is not None): fpdc1 = copy(fpdc1) + '.bkp'
                        if(fpdc2 is not None): fpdc2 = copy(fpdc2) + '.bkp'
                        if(fpA is not None): fpA = copy(fpA) + '.bkp'
                        
                        bkpFileToBeSaved = False
                    else: bkpFileToBeSaved = True
                    
                if(((_char == 6) or (_char == 61) or (_char == 7)) and (_n_mpi <= 1)): #OC20062021 (copy / update CSD only if total distribution is required)
                #if((_char == 6) and (_n_mpi <= 1)): #OC03032021 (copy / update CSD only if total distribution is required)
                #if(_char == 6): #OC18022021
                    srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [4]) #Filling-in "symmetrical" part of the Hermitian Mutual Intensity distribution
                    
                #if(_char == 40): #OC03052018
                if((_char == 40) or (_char == 41)): #OC13072019
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0)
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, fp, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0) #OC14082018
                    srwl_uti_save_intens(resStokes.arS, meshRes, fp, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0, _form = _file_form) #OC17072021
                    if(resStokes2 is not None):
                        #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) 
                        #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC060502018
                        if(_char == 40): #OC13072019
                            srwl_uti_save_intens(resStokes2.arS, resStokes2.mesh, fp1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1, _form = _file_form) #OC17072021
                            #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, fp1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC14082018
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC14082018
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs X

                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs X
                        #OCTEST14112020
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(_rel_zer_tol=0), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs X
                        srwl_uti_save_intens(resStokes2.to_deg_coh(_rel_zer_tol=0), resStokes2.mesh, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs X

                    if(resStokes3 is not None):
                        #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) 
                        #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC060502018
                        if(_char == 40): #OC13072019
                            srwl_uti_save_intens(resStokes3.arS, resStokes3.mesh, fp2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1, _form = _file_form) #OC17072021
                            #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, fp2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC14082018
                        #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                        #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC14082018
                        #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs Y

                        #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs Y
                        #OCTEST14112020
                        #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(_rel_zer_tol=0), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs Y
                        srwl_uti_save_intens(resStokes3.to_deg_coh(_rel_zer_tol=0), resStokes3.mesh, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs Y

                #elif(_char == 4): #OC03052018
                elif((_char == 4) or (_char == 5)): #OC13072019
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC060502018
                    if(_char == 4): #OC13072019
                        srwl_uti_save_intens(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                    #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                    #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC14082018
                    #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs X
                    #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs X
                    srwl_uti_save_intens(resStokes.to_deg_coh(), meshRes, fpdc1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs X
                    if((resStokes2 is not None) and (meshRes2 is not None)):
                        #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) 
                        #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC060502018
                        if(_char == 4): #OC13072019
                            srwl_uti_save_intens(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021
                            #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0)
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0) #OC14082018
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = doMutual, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs Y
                        #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs Y
                        srwl_uti_save_intens(resStokes2.to_deg_coh(), meshRes2, fpdc2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs Y
                else:
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC30052017
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                    #srwl_uti_save_intens(resStokes.arS, meshRes, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC18022021
                    #srwl_uti_save_intens(resStokes.arS, resStokes.mesh, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17042021
                    srwl_uti_save_intens(resStokes.arS, resStokes.mesh, fp1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form, _wfr = wfrA) #OC18062021
                    #To use the above function everywhere

                    if((resStokes2 is not None) and (meshRes2 is not None)): #OC30052017
                        #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) 
                        #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                        srwl_uti_save_intens(resStokes2.arS, meshRes2, fp2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021

                if(_pres_ang == 2): #OC24122018
                    srwl_uti_save_intens(resStokesA.arS, meshResA, fpA, numComp, _arLabels = resLabelsToSaveA, _arUnits = resUnitsToSaveA, _mutual = 0, _cmplx = 0, _form = _file_form) #OC17072021
                    #srwl_uti_save_intens_ascii(resStokesA.arS, meshResA, fpA, numComp, _arLabels = resLabelsToSaveA, _arUnits = resUnitsToSaveA, _mutual = 0, _cmplx = 0)

                #DEBUG
                #srwl_uti_save_intens_ascii(workStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017
                #END DEBUG
                #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                
                iSave = 0
    #DEBUG
    #srwl_uti_save_text("Exiting srwl_wfr_emit_prop_multi_e", _file_path + "." + str(rank) + "e.dbg")
    #END DEBUG

    #OC04042021
    #DEBUG (Barrier commented-out)
    if(nProc > 1): 
        if(rank != rankMaster): #OC26042021
            comMPI.Barrier() #Attempt to prevent crashes because some processes quit too early

    if((rank == rankMaster) or (nProc == 1)): #OC02032021
    #if((rank == 0) or (nProc == 1)):
        #Saving final results:
        #if(_file_path != None):
        #if(file_path1 is not None): #OC30052017

        #print('srwl_uti_save_intens_ascii ... ', end='') #DEBUG
        #t0 = time.time(); #DEBUG
            
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC26042016
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017

        if((nProc == 1) or ((rank == rankMaster) and (iSave > 0))): #OC22042021 (save only if necessary)
        #if((nProc == 1) or ((rank == rankMaster) and (iSave > 0)) or ((rank == 0) and (_char == 6) and (_n_mpi > 1))): #OC21042021 (save only if necessary)

            if(nProc > 1): #OC31102021
                particle_number = (i + 1)*_n_part_avg_proc
                srwl_uti_save_stat_wfr_emit_prop_multi_e(particle_number, actNumPartTot, filename=log_path) #OC31102021

            #if(_char == 40): #OC03052018
            if((_char == 40) or (_char == 41)): #OC13072019
                srwl_uti_save_intens(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0, _form = _file_form) #OC17072021
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0)
                if(resStokes2 is not None):
                    #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                    if(_char == 40): #OC13072019
                        srwl_uti_save_intens(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1, _form = _file_form) #OC17072021
                        #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs X

                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs X
                    #OCTEST14112020
                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(_rel_zer_tol=0), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs X
                    srwl_uti_save_intens(resStokes2.to_deg_coh(_rel_zer_tol=0), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs X

                if(resStokes3 is not None):
                    #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                    if(_char == 40): #OC13072019
                        srwl_uti_save_intens(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1, _form = _file_form) #OC17072021
                        #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                    #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                    #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs Y

                    #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs Y
                    #OCTEST14112020
                    #srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(_rel_zer_tol=0), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs Y
                    srwl_uti_save_intens(resStokes3.to_deg_coh(_rel_zer_tol=0), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs Y

            #elif(_char == 4): #OC03052018
            elif((_char == 4) or (_char == 5)): #OC13072019
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                if(_char == 4): #OC13072019
                    srwl_uti_save_intens(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 1, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs X
                #srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs X
                srwl_uti_save_intens(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs X
                if((resStokes2 is not None) and (meshRes2 is not None)):
                    #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                    if(_char == 4): #OC13072019
                        srwl_uti_save_intens(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021
                        #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0)
                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = doMutual, _cmplx = 0) #OC12072019 # Deg. of Coh. Cut vs Y
                    #srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0) #OC16072019 # Deg. of Coh. Cut vs Y
                    srwl_uti_save_intens(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSaveDC, _arUnits = resUnitsToSaveDC, _mutual = 2, _cmplx = 0, _form = _file_form) #OC17072021 # Deg. of Coh. Cut vs Y
            else:
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC16012017

                if((nProc == 1) and ((_char == 6) or (_char == 61) or (_char == 7))): #OC20062021
                #if((nProc == 1) and (_char == 6)): #OC18062021
                    srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [4]) #Fill-in "symmetrical" part of the Hermitian MI distribution again (before final saving)

                    if((_char == 61) or (_char == 7)): #OC27062021
                        cohModes, eigVals = srwl_wfr_cmd(resStokes, _n_modes=_n_cm, _awfr=wfrA)
                        #DEBUG OC30042022
                        #cohModes, eigVals = srwl_wfr_cmd(resStokes, _n_modes=_n_cm, _awfr=None)
                        #END DEBUG

                        #DEBUG
                        t0 = time.time()
                        #END DEBUG

                        #file_path_cm = srwl_wfr_fn(_file_path, _type=7, _form='hdf5')
                        #srwl_uti_save_wfr_cm_hdf5(cohModes, _awfr=wfrA, _file_path=fp_cm)
                        srwl_uti_save_wfr_cm_hdf5(cohModes, None, _awfr=wfrA, _file_path=fp_cm) #OC28062021

                        #DEBUG
                        t1 = round(time.time() - t0, 3)
                        print('Coherent Modes file was saved in:', t1, 's')
                        sys.stdout.flush()
                        #END DEBUG

                        if(_char == 7): #OC02072021: delete CSD / MI file(?)
                            if os.path.exists(file_path1): os.remove(file_path1)
                            return cohModes, eigVals, resStokes.mesh 
                        else:
                            return resStokes, cohModes, eigVals #?

                srwl_uti_save_intens(resStokes.arS, resStokes.mesh, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form, _wfr = wfrA) #OC18062021
                #srwl_uti_save_intens(resStokes.arS, resStokes.mesh, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC21042021
                #srwl_uti_save_intens(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC18022021
                #if(_file_form == 'ascii'): #OC05022021
                #    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC14082018
                #elif(_file_form == 'hdf5'): #OC05022021
                #    srwl_uti_save_intens_hdf5(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))

                if((resStokes2 is not None) and (meshRes2 is not None)): #OC03052018
                    srwl_uti_save_intens(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC17072021
                    #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC16012017

            if(_pres_ang == 2): #OC24122018
                srwl_uti_save_intens(resStokesA.arS, meshResA, file_pathA, numComp, _arLabels = resLabelsToSaveA, _arUnits = resUnitsToSaveA, _mutual = 0, _cmplx = 0, _form = _file_form) #OC17072021
                #srwl_uti_save_intens_ascii(resStokesA.arS, meshResA, file_pathA, numComp, _arLabels = resLabelsToSaveA, _arUnits = resUnitsToSaveA, _mutual = 0, _cmplx = 0)

        if(nProc > 1): #OC26042021
            comMPI.Barrier() #Attempt to prevent crashes because some processes quit too early

            #OC18062021 (do this only if nProc > 1)
            if(((_char == 6) or (_char == 61) or (_char == 7)) and (rank == 0)): #OC28102021
            #if(((_char == 6) or (_char == 61) or (_char == 7)) and (rank == 0) and (_n_mpi > 1)): #OC20062021
            #if((_char == 6) and (rank == 0) and (_n_mpi > 1)): #OC22042021 moved here from top (assemble the final total CSD from different parts stored in files, by the Master of the first group (rank = 0))
                #Summing-up/averaging partial CSDs

                ##del resStokes
                #resStokes0 = resStokes #OC23042021
                #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual, _n_comp = 1) #OC04032021 (to store final MI)
                #OC25042021: Commented the above out for the version/case when total CSD is calculated by each CSD group, and the final CSD is averaged

                if(_n_mpi > 1): #OC28102021
                    doFinSave = False if(_char == 7) else True
                    CSD, wfrA_dummy = srwl_wfr_csd_avg(fpCore, fpExt, 0, _n_mpi-1, _csd0=resStokes, _awfr=wfrA, _form=_file_form, _do_fin_save=doFinSave, _do_del_aux_files=_del_aux_files, _do_sum=doPropCM) #OC31102021
                    #CSD, wfrA_dummy = srwl_wfr_csd_avg(fpCore, fpExt, 0, _n_mpi-1, _csd0=resStokes, _awfr=wfrA, _form=_file_form, _do_fin_save=doFinSave, _do_del_aux_files=_del_aux_files) #OC10102021
                else:
                    CSD = resStokes #OC28102021

                #CSD = srwl_wfr_csd_avg(fpCore, fpExt, 0, _n_mpi-1, _csd0=resStokes, _awfr=wfrA, _form=_file_form, _do_fin_save=doFinSave, _do_del_aux_files=_del_aux_files) #OC27062021

                # #DEBUG
                # print('Rank=', rank, ' About to start collecting CSD data from different MPI Groups')
                # sys.stdout.flush()
                # #END DEBUG

                # for j in range(1, _n_mpi): #OC25042021
                # #for j in range(_n_mpi):

                #     #if(j > 0): #OC25042021 (commented-out)
                #     curFP = fpCore + '_' + repr(j) + fpExt
                #     curFP = srwl_wfr_fn(curFP, 3) #adds suffix "_mi" before extension

                #     resPartMI = srwl_uti_read_intens(_file_path=curFP, _form=_file_form)
                #     arPartMI = resPartMI[0]
                #     meshPartMI = resPartMI[1]

                #     #else: #OC25042021 (commented-out)
                #     #    arPartMI = resStokes0.arS
                #     #    meshPartMI = resStokes0.mesh

                #     #DEBUG
                #     print('Rank=', rank, ' Adding CSD data from MPI Group:', j, ' is about to start')
                #     sys.stdout.flush()
                #     #END DEBUG

                #     #DEBUG
                #     #print('Rank=', rank, ' Some resStokes.arS data BEFORE averaging with data from MPI Group:', j)
                #     #print('resStokes.arS[0]=', resStokes.arS[0])
                #     #print('resStokes.arS[2*nx*ny + 1]=', resStokes.arS[2*resStokes.mesh.nx*resStokes.mesh.ny + 1])
                #     #print('resStokes.arS[(2*nx*ny)*2 + 3]=', resStokes.arS[(2*resStokes.mesh.nx*resStokes.mesh.ny)*2 + 3])
                #     #END DEBUG

                #     #DEBUG
                #     #print('Rank=', rank, ' Some arPartMI data BEFORE averaging it with resStokes.arS data, j=', j)
                #     #print('arPartMI[0]=', arPartMI[0])
                #     #print('arPartMI[2*nx*ny + 1]=', arPartMI[2*meshPartMI.nx*meshPartMI.ny + 1])
                #     #print('arPartMI[(2*nx*ny)*2 + 3]=', arPartMI[(2*meshPartMI.nx*meshPartMI.ny)*2 + 3])
                #     #END DEBUG

                #     #OC25042021
                #     srwl.UtiIntProc(resStokes.arS, resStokes.mesh, arPartMI, meshPartMI, [1, j]) #Average new total MI (arPartMI) with the total MI (resStokes.arS), using information in meshPartMI
                #     #srwl.UtiIntProc(resStokes.arS, resStokes.mesh, arPartMI, meshPartMI, [1]) #Add partial MI (arPartMI) to the total MI (resStokes.arS), using information in meshPartMI

                #     #DEBUG
                #     print('Rank=', rank, ' CSD data from MPI Group:', j, ' was added')
                #     sys.stdout.flush()
                #     #END DEBUG

                #     #DEBUG
                #     #print('Rank=', rank, ' Some resStokes.arS data AFTER averaging with data from MPI Group:', j)
                #     #print('resStokes.arS[0]=', resStokes.arS[0])
                #     #print('resStokes.arS[2*nx*ny + 1]=', resStokes.arS[2*resStokes.mesh.nx*resStokes.mesh.ny + 1])
                #     #print('resStokes.arS[(2*nx*ny)*2 + 3]=', resStokes.arS[(2*resStokes.mesh.nx*resStokes.mesh.ny)*2 + 3])
                #     #END DEBUG

                #     #OC25042021: Commented-out the line below, since all CSDs are avaraged with the one of the first MPI group
                #     #if(j == 0): del resStokes0

                # srwl.UtiIntProc(resStokes.arS, resStokes.mesh, None, None, [4]) #Fill-in "symmetrical" part of the Hermitian MI distribution

                # #DEBUG
                # print('Rank=', rank, ' Symmetrical parts of the CSD data (Hermitian matrix) were filled-out')
                # sys.stdout.flush()
                # #END DEBUG

                # #Final saving is done here
                # #OC22042021
                # fpResMI = srwl_wfr_fn(fpCore + fpExt, 3) #adds suffix "_mi" before extension

                # srwl_uti_save_intens(resStokes.arS, resStokes.mesh, fpResMI, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form, _wfr = wfrA) #OC18062021
                # #srwl_uti_save_intens(resStokes.arS, resStokes.mesh, fpResMI, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0), _form = _file_form) #OC21042021

                # #DEBUG
                # print('Rank=', rank, ' Total CSD data saved')
                # sys.stdout.flush()
                # #END DEBUG

                if((_char == 61) or (_char == 7)): #OC27062021
                    cohModes, eigVals = srwl_wfr_cmd(CSD, _n_modes=_n_cm, _awfr=wfrA)

                    #DEBUG
                    t0 = time.time()
                    #END DEBUG

                    #file_path_cm = srwl_wfr_fn(_file_path, _type=7, _form='hdf5')
                    #srwl_uti_save_wfr_cm_hdf5(cohModes, _awfr=wfrA, _file_path=fp_cm)
                    srwl_uti_save_wfr_cm_hdf5(cohModes, None, _awfr=wfrA, _file_path=fp_cm) #OC28062021

                    #DEBUG
                    t1 = round(time.time() - t0, 3)
                    print('Coherent Modes file was saved in:', t1, 's')
                    sys.stdout.flush()
                    #END DEBUG

                    if(_char == 7): #OC02072021: delete CSD / MI file(?)
                        if os.path.exists(file_path1): os.remove(file_path1)
                        return cohModes, eigVals, CSD.mesh
                    else: 
                        return CSD, cohModes, eigVals #?

        #DEBUG
        #print('Rank=', rank, ' Reached Barrier before Return')
        #sys.stdout.flush()
        #END DEBUG

        #OC22042021 (?)
        #OC24042021 (commented-out)
        #if(nProc > 1): comMPI.Barrier() #Attempt to prevent crashes because some processes quit too early

        #DEBUG
        #print('Rank=', rank, ' PASSED Barrier, about to Return')
        #sys.stdout.flush()
        #END DEBUG

        if(rank != 0): return None #OC23042021

        #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
        #return resStokes
        if(resStokes2 is None): #OC30052017
            #return resStokes
            if(resStokesA is None): return resStokes #OC24122018
            else: return resStokes, resStokesA
        elif(resStokes3 is None): #OC03052018
            #return resStokes, resStokes2
            if(resStokesA is None): return resStokes, resStokes2
            else: return resStokes, resStokes2, resStokesA
        else:
            #return resStokes, resStokes2, resStokes3
            if(resStokesA is None): return resStokes, resStokes2, resStokes3
            else: return resStokes, resStokes2, resStokes3, resStokesA
    else:

        #DEBUG
        #print('Rank=', rank, ' Reached Barrier before Return')
        #sys.stdout.flush()
        #END DEBUG

        #OC24042021 (commented-out)
        #if(nProc > 1): comMPI.Barrier() #Attempt to prevent crashes because some processes quit too early

        #DEBUG
        #print('Rank=', rank, ' PASSED Barrier, about to Return')
        #sys.stdout.flush()
        #END DEBUG

        return None

#****************************************************************************
#****************************************************************************
#Import of modules requiring classes defined in this smodule
#****************************************************************************
#****************************************************************************
try: #OC04032024
    from .srwl_uti_src import *
except:
    from srwl_uti_src import *

#from srwl_uti_src import *

#****************************************************************************
#****************************************************************************
# Help to main functions implemented in C/C++ (available through srwlpy.pyd/.so)
#****************************************************************************
#****************************************************************************
helpCalcMagnField = """CalcMagnField(_outMagFld3DC, _inMagFldC, _inPrec)
function calculates (tabulates) 3D magnetic field created by different magnetic field sources / elements
:param _outMagFld3DC: output magnetic field container (instance of SRWLMagFldC) with the tabulated 3D magnetic field element
       (instance of SRWLMagFld3D) 
:param _inMagFldC: input magnetic field container (instance of SRWLMagFldC) of magnetic field sources / elements
:param _inPrec: optional array of precision parameters
        _precPar[0] defines the type of calculation: =0 -standard calculation, =1 -interpolation vs one parameter, =2 -interpolation vs two parameters
        _precPar[1]: first parameter value the field has to be interpolated for
        _precPar[2]: second parameter value the field has to be interpolated for
        _precPar[3]: specifies type of interpolation: =1 -(bi-)linear, =2 -(bi-)quadratic, =3 -(bi-)cubic 
"""
helpCalcPartTraj = """CalcPartTraj(_prtTrj, _inMagFldC, _inPrec)
function calculates charged particle trajectory in external 3D magnetic field (in Cartesian laboratory frame)
:param _prtTrj: input / output trajectory structure (instance of SRWLPrtTrj);
       note that all data arrays should be allocated in Python script before calling this function;
       initial conditions and particle type must be specified in _prtTrj.partInitCond;
       the initial conditions are assumed to be given for ct = 0,
       however the trajectory will be calculated for the mesh defined by _prtTrj.np, _prtTrj.ctStart, _prtTrj.ctEnd
:param _inMagFldC: input magnetic field container structure (instance of SRWLMagFldC)
:param _inPrec: input list of calculation method ID and precision parameters;
       _inPrec[0]: integration method ID:
                   =1 -use the fourth-order Runge-Kutta (R-K), with the precision driven by number of points
                   =2 -use the fifth-order R-K
       _inPrec[1],[2],[3],[4],[5]: optional absolute precision values for X[m],X'[rad],Y[m],Y'[rad],Z[m] 
                   to be taken into account only for R-K fifth order or higher (yet to be tested!!)
       _inPrec[6]: tolerance (default = 1) for R-K fifth order or higher
       _inPrec[7]: maximal number of auto-steps for R-K fifth order or higher (default = 5000)
"""
helpCalcPartTrajFromKickMatr = """CalcPartTrajFromKickMatr(_prtTrj, _inKickM, _inPrec)
function calculates charged particle trajectory from one or a list of kick-matrices
:param _prtTrj: input / output trajectory structure (instance of SRWLPrtTrj);
       note that all data arrays should be allocated in Python script before calling this function;
       initial conditions and particle type must be specified in _partTraj.partInitCond;
       the initial conditions are assumed to be given for ct = 0,
       however the trajectory will be calculated for the mesh defined by _prtTrj.np, _prtTrj.ctStart, _prtTrj.ctEnd
:param _inKickM: input kick-matrix (instance of SRWLKickM) or a list of such kick-matrices
:param _inPrec: input list of calculation parameters:
       _inPrec[0]: switch specifying whether the new trajectory data should be added to pre-existing trajectory data (=1, default)
       or it should override any pre-existing trajectory data (=0)
"""
helpCalcElecFieldSR = """CalcElecFieldSR(_wfr, _inPrtTrj, _inMagFldC, _inPrec)
function calculates Electric Field (Wavefront) of Synchrotron Radiation by a relativistic charged particle
traveling in external 3D magnetic field
:param _wfr: input / output resulting Wavefront structure (instance of SRWLWfr);
       all data arrays should be allocated in Python script before calling this function;
       the emitting particle beam, radiation mesh, presentation, etc., should be specified in this structure at input
:param _inPrtTrj: optional input pre-calculated particle trajectory structure (instance of SRWLPrtTrj);
       the initial conditions and particle type must be specified in _inPrtTrj.partInitCond;
       if the trajectory data arrays (_inPrtTrj.arX, _inPrtTrj.arXp, _inPrtTrj.arY, _inPrtTrj.arYp) are defined,
       the SR will be calculated from these data; if these arrays are not defined, or if _inPrtTrj =0, the function will attempt
       to calculate the SR from the magnetic field data (_inMagFldC) which has to be supplied
:param _inMagFldC: optional input magnetic field (container) structure (instance of SRWLMagFldC);
       to be taken into account only if particle trajectroy arrays (_inPrtTrj.arX, _inPrtTrj.arXp, _inPrtTrj.arY, _inPrtTrj.arYp)
       are not defined
:param _inPrec: input list of precision parameters:
       _inPrec[0]: method ID: =0 -"manual", =1 -"auto-undulator", =2 -"auto-wiggler")
       _inPrec[1]: step size (for "manual" method, i.e. if _inPrec[0]=0) or relative precision
                   (for "auto-undulator" or "auto-wiggler" methods, i.e. if _inPrec[0]=1 or _inPrec[0]=2)
       _inPrec[2]: longitudinal position [m] to start integration (effective if _inPrec[2] < _inPrec[3])
       _inPrec[3]: longitudinal position [m] to finish integration (effective if _inPrec[2] < _inPrec[3])
       _inPrec[4]: number of points to use for trajectory calculation 
       _inPrec[5]: calculate terminating terms or not:
                   =0 -don't calculate two terms,
                   =1 -do calculate two terms,
                   =2 -calculate only upstream term,
                   =3 -calculate only downstream term
       _inPrec[6]: sampling factor (for propagation, effective if > 0)
"""
helpCalcElecFieldGaussian = """CalcElecFieldGaussian(_wfr, _inGsnBm, _inPrec)
function calculates Electric Field (Wavefront) of a coherent Gaussian beam
:param _wfr: input / output resulting Wavefront structure (instance of SRWLWfr);
       all data arrays should be allocated in Python script before calling this function;
       the emitting particle beam, radiation mesh, presentation, etc., should be specified in this structure at input
:param _inGsnBm: input coherent Gaussian beam parameters structure (instance of SRWLGsnBm)
:param _inPrec: input list of precision parameters:
       _inPrec[0]: sampling factor (for propagation, effective if > 0)
"""
helpCalcElecFieldPointSrc = """CalcElecFieldPointSrc(_wfr, _inPtSrc, _inPrec)
function calculates Electric Field (Wavefront) of a Pont Source (i.e. spherical wave)
:param _wfr: input / output resulting Wavefront structure (instance of SRWLWfr); 
        all data arrays should be allocated in a calling function/application; 
        the mesh, presentation, etc., should be specified in this structure at input
:param _inPtSrc: input Point Source parameters structure (instance of SRWLPtSrc)
:param _inPrec: input list of precision parameters:
       _inPrec[0]: sampling factor (for propagation, effective if > 0)
"""
helpCalcStokesUR = """CalcStokesUR(_stk, _inElBeam, _inUnd, _inPrec)
function calculates Stokes parameters of Undulator Radiation (UR) by a relativistic finite-emittance electron beam
traveling in periodic magnetic field of an undulator
:param _stk: input / output resulting Stokes structure (instance of SRWLStokes);
       all data arrays should be allocated in Python script before calling this function; the mesh, presentation, etc.,
       should be specified in this structure at input
:param _inElBeam: input electron beam structure (instance of SRWLPartBeam)
:param _inUnd: input undulator (periodic magnetic field) structure (instance of SRWLMagFldU)
:param _inPrec: input list of precision parameters:
       _inPrec[0]: initial harmonic of UR spectrum
       _inPrec[1]: final harmonic of UR spectrum
       _inPrec[2]: longitudinal integration precision parameter (nominal value is 1.0, for better accuracy make it > 1.0)
       _inPrec[3]: azimuthal integration precision parameter (nominal value is 1.0, for better accuracy make it > 1.0)
       _inPrec[4]: calculate flux (=1) or intensity (=2)
"""
helpCalcPowDenSR = """CalcPowDenSR(_stk, _inElBeam, _inPrtTrj, _inMagFldC, _inPrec)
function calculates Power Density distribution of Synchrotron Radiation by a relativistic finite-emittance electron beam
traveling in arbitrary magnetic field
:param _stk: input / output resulting Stokes structure (instance of SRWLStokes); 
       all data arrays should be allocated in Python script before calling this function; the mesh, presentation, etc.,
       should be specified in this structure at input; the Power Density data will be written to _stk.arS
:param _inElBeam: input electron beam structure (instance of SRWLPartBeam)
:param _inPrtTrj: input trajectory structure (instance of SRWLPrtTrj);
       can be =0; in such case, the power density is calculated based on _inElBeam and _inMagFldC
:param _inMagFldC: input magnetic field container structure (instance of SRWLMagFldC);
       can be =0; in such case, power density is calculated from _inPrtTrj (if _inPrtTrj != 0) and _inElBeam (if _inElBeam != 0))
:param _inPrec: input list of precision parameters:
       _inPrec[0]: precision factor (=1.0 default, >1.0 for more precision)
       _inPrec[1]: power density computation method (=1 -"near field" (default), =2 -"far field")
       _inPrec[2]: initial longitudinal position [m] (effective if < _inPrec[3])
       _inPrec[3]: final longitudinal position [m] (effective if > _inPrec[2])
       _inPrec[4]: number of points to use for trajectory calculation 
"""
helpCalcIntFromElecField = """CalcIntFromElecField(_arI, _inWfr, _inPol, _inIntType, _inDepType, _inE, _inX, _inY)
function calculates/"extracts" Intensity from pre-calculated Electric Field
:param _arI: output resulting Intensity array (should be allocated in Python script before calling this function)
:param _inWfr: input pre-calculated Wavefront structure (instance of SRWLWfr)
:param _inPol: input switch specifying polarization component to be extracted:
               =0 -Linear Horizontal; 
               =1 -Linear Vertical; 
               =2 -Linear 45 degrees; 
               =3 -Linear 135 degrees; 
               =4 -Circular Right; 
               =5 -Circular Left; 
               =6 -Total
:param _inIntType: input switch specifying "type" of a characteristic to be extracted:
               =0 -"Single-Electron" Intensity; 
               =1 -"Multi-Electron" Intensity; 
               =2 -"Single-Electron" Flux; 
               =3 -"Multi-Electron" Flux; 
               =4 -"Single-Electron" Radiation Phase; 
               =5 -Re(E): Real part of Single-Electron Electric Field;
               =6 -Im(E): Imaginary part of Single-Electron Electric Field;
               =7 -"Single-Electron" Intensity, integrated over Time or Photon Energy (i.e. Fluence)
               =8 -"Single-Electron" Mutual Intensity (i.e. E(r)E*(r'))
               under development: =9 -"Multi-Electron" Mutual Intensity 
:param _inDepType: input switch specifying type of dependence to be extracted:
               =0 -vs e (photon energy or time);
               =1 -vs x (horizontal position or angle);
               =2 -vs y (vertical position or angle);
               =3 -vs x&y (horizontal and vertical positions or angles);
               =4 -vs e&x (photon energy or time and horizontal position or angle);
               =5 -vs e&y (photon energy or time and vertical position or angle);
               =6 -vs e&x&y (photon energy or time, horizontal and vertical positions or angles);
:param _inE: input photon energy [eV] or time [s] to keep fixed (to be taken into account for dependences vs x, y, x&y)
:param _inX: input horizontal position [m] to keep fixed (to be taken into account for dependences vs e, y, e&y)
:param _inY: input vertical position [m] to keep fixed (to be taken into account for dependences vs e, x, e&x)
"""
helpCalcTransm = """CalcTransm(_opT, _inDelta, _inAttenLen, _inObjShapeDefs, _inPrec)
Sets Up Transmittance for an Optical Element defined from a list of 3D (nano-) objects, e.g. for simulating samples for coherent scattering experiments
:param _opT: input/output Optical Transmission object to set up
:param _inDelta: input array of (spectral) Refractive Index Decrement data
:param _inAttenLen input array of (spectral) Attenuation Length data
:param _inObjShapeDefs input list of 3D object shape definitions
:param _inPrec input array of precision parameters (currently unused)
"""
helpResizeElecField = """ResizeElecField(_wfr, _inType, _inPar)
function resizes Electric Field Wavefront vs transverse positions / angles or photon energy / time
:param _wfr: input / output Wavefront structure (instance of SRWLWfr)
:param _inType: input character specifying whether the resizing should be done vs positions / angles ('c')
       or vs photon energy / time ('f')
:param _inPar: input list of parameters (meaning depends on value of _inType):
       if(_inType == 'c'):
           _inPar[0]: method (=0 -regular method, without FFT, =1 -"special" method involving FFT)
           _inPar[1]: range resizing factor for horizontal position / angle
                      (range will be decreased if 0 < _inPar[1] < 1. and increased if _inPar[1] > 1.)
           _inPar[2]: resolution resizing factor for horizontal position / angle
                      (resolution will be decreased if 0 < _inPar[2] < 1. and increased if _inPar[2] > 1.)
           _inPar[3]: range resizing factor for vertical position / angle
                      (range will be decreased if 0 < _inPar[3] < 1. and increased if _inPar[3] > 1.)
           _inPar[4]: resolution resizing factor for vertical position / angle
                      (resolution will be decreased if 0 < _inPar[4] < 1. and increased if _inPar[4] > 1.)
           _inPar[5]: relative horizontal wavefront center position / angle at resizing
                      (default is 0.5; in that case the resizing will be symmetric)
           _inPar[6]: relative vertical wavefront center position / angle at resizing
                      (default is 0.5; in that case the resizing will be symmetric)
       if(_inType == 'f'):
           _inPar[0]: method (=0 -regular method, without FFT, =1 -"special" method involving FFT)
           _inPar[1]: range resizing factor for photon energy / time
                      (range will be decreased if 0 < _inPar[1] < 1. and increased if _inPar[1] > 1.)
           _inPar[2]: resolution resizing factor for photon energy / time
                      (resolution will be decreased if 0 < _inPar[2] < 1. and increased if _inPar[2] > 1.)
           _inPar[3]: relative photon energy / time center position at resizing
                      (default is 0.5; in that case the resizing will be symmetric)
"""
helpSetRepresElecField = """SetRepresElecField(_wfr, _inRepr)
function changes Representation of Electric Field: positions <-> angles, frequency <-> time
:param _wfr: input / output Wavefront structure (instance of SRWLWfr)
:param _inRepr: input character specifying desired representation:
                ='c' for coordinate, ='a' for angle, ='f' for frequency, ='t' for time
"""
helpPropagElecField = """PropagElecField(_wfr, _inOptC)
function propagates Electric Field Wavefront through Optical Elements and free space
:param _wfr: input / output Wavefront structure (instance of SRWLWfr)
:param _inOptC: input container of optical elements (instance of SRWLOptC) the propagation should be done through;
       note that lists of optical elements and the corresponding propagation parameters have to be defined in
       _inOptC.arOpt and _inOptC.arProp respectively (see help/comments to SRWLOptC class)
"""
helpUtiFFT = """UtiFFT(_data, _mesh, _inDir)
function performs 1D or 2D in-place Fast Fourier Transform (as defined by arguments)
:param _data: input / output float (single-precision) type array of data to be Fourier-transformed;
       in the case of a 2D transform, the data should be "C-aligned" in 1D array,
       with the first dimension being "inner" (i.e. most frequently changing);
       the FFT is performed "in place", i.e. the input _data will be replaced by a resulting data
:param _mesh: input / output list specifying (equidistant, regular) mesh of the data to be transformed:
       _mesh[0]: start value of the first argument
       _mesh[1]: step size value of the first argument
       _mesh[2]: number of points over the first argument (note: it should be even!)
       _mesh[3]: (optional, to be used for 2D FFT) start value of the second argument
       _mesh[4]: (optional, to be used for 2D FFT) step size value of the second argument
       _mesh[5]: (optional, to be used for 2D FFT) number of points of the second argument (note: it should be even!)
       if len(_mesh) == 3, 1D FFT will be performed
       else if len(_mesh) == 6, 2D FFT will be performed
       the input _mesh will be replaced by a resulting mesh
:param _inDir: input integer number specifying FFT "direction": >0 means forward FFT, <0 means backward FFT
"""
helpUtiConvWithGaussian = """UtiConvWithGaussian(_data, _inMesh, _inSig)
function performs convolution of 1D or 2D data wave with 1D or 2D Gaussian (as defined by arguments)
:param _data: input / output float (single-precision) type array of data to be convolved with Gaussian;
       in the case of a 2D convolution, the data should be "C-aligned" in 1D array,
       with the first dimension being "inner" (i.e. most frequently changing);
       the convolution is performed "in place", i.e. the input _data will be replaced by a resulting data
:param _inMesh: input list specifying (equidistant, regular) mesh of the data to be transformed:
       _inMesh[0]: start value of the first argument
       _inMesh[1]: step size value of the first argument
       _inMesh[2]: number of points over the first argument
       _inMesh[3]: (optional, to be used for 2D convolution) start value of the second argument
       _inMesh[4]: (optional, to be used for 2D convolution) step size value of the second argument
       _inMesh[5]: (optional, to be used for 2D convolution) number of points of the second argument
       if len(_mesh) == 3, 1D convolution will be performed
       else if len(_mesh) == 6, 2D convolution will be performed
:param _inSig: input list of central 2nd order statistical moments of 1D or 2D Gaussian,
       and possibly a coefficient before cross-term:
       _inSig[0]: RMS size of the Gaussian in first dimension
       _inSig[1]: (optional) RMS size of a 2D Gaussian in second dimension
       _inSig[2]: (optional) coefficient before cross-term in exponent argument of a 2D Gaussian
       i.e. _inSig[] = [sigX, sigY, alp] defines a "tilted" normalized 2D Gaussian (vs x, y): 
       (sqrt(1 - (alp*sigX*sigY)**2)/(2*Pi*sigX*sigY))*exp(-x**2/(2*sigX**2) - y**2/(2*sigY^2) - alp*x*y)
"""
helpUtiIntInf = """UtiIntInf(_inData, _inMesh, _inPrec)
function calculates basic statistical characteristics of an intensity distribution
:param _inData flat (C-aligned) array of intensity distribution data to be analyzed
:param _inMesh instance of SRWLRadMesh describing mesh (grid) of radiation intensity distribution
:param _inPrec optional array / list of precision parameters:
       _inPrec[0]: method to be used for determining FWHM values: =0 means start intensity scan from extremities of the distribution, =1 means start intensity scan from the peak of the distribution
       _inPrec[1]: (optional) fraction of maximum for determining additional Full Width at a Fraction of Maximum value over 1st dimension
       _inPrec[2]: (optional) fraction of maximum for determining additional Full Width at a Fraction of Maximum value over 2nd dimension
       _inPrec[3]: (optional) fraction of maximum for determining additional Full Width at a Fraction of Maximum value over 3rd dimension
:return list of distribution characteristics to be calculated:
       [0]: peak (max.) intensity
       [1]: position of peak intensity vs 1st dimension
       [2]: position of peak intensity vs 2nd dimension
       [3]: position of peak intensity vs 3rd dimension (reserved for future use)
       [4]: FWHM value of intensity distribution vs 1st dimension
       [5]: FWHM value of intensity distribution vs 2nd dimension
       [6]: FWHM value of intensity distribution vs 3rd dimension (reserved for future use)
       [7]: (optional) additional Full Width at a Fraction of Maximum value of intensity distribution over 1st dimension (the fraction is defined by _inPrec[1])
       [8]: (optional) additional Full Width at a Fraction of Maximum value of intensity distribution over 2nd dimension (the fraction is defined by _inPrec[2])
       [9]: (optional) additional Full Width at a Fraction of Maximum value of intensity distribution over 3rd dimension (the fraction is defined by _inPrec[3])
       more characteristics to be added
"""
helpUtiIntProc = """UtiIntProc(_data, _mesh, _inData, _inMesh, _inPrec)
function performs misc. operations on intensity distribution (or similar C-aligned) arrays
:param _data input / output flat (C-aligned) array of intensity distribution data
:param _mesh input / output instance of SRWLRadMesh describing mesh (grid) of radiation intensity distribution _data
:param _inData input flat (C-aligned) array of intensity distribution data to be processed
:param _inMesh input instance of SRWLRadMesh describing mesh (grid) of radiation intensity distribution _inData to be processed
:param _inPrec array / list of precision parameters:
       _inPrec[0]: defines type of the operation and the meaning of other elements dependent on it:
               =1 -add (without or with averaging) intensity or mutual intensity distribution _inData to the distribution _data and store result in _data (depending on the meshes of the two distributions, _inMesh and _mesh, it may or may not do interpolation of _inData)
                   in that case, the meaning of the subsequent parameters stored in _inPrec is:
                   _inPrec[1] defines whether simple summation (=-1, default), or averaging should take place, in the latter case it is the iteration number (>=0)
               =2 -find average of intensity distribution _inData and the distribution _data, assuming it to be a given iteration, and store result in _data (depending on the meshes of the two distributions, _inMesh and _mesh, it may or mey not do interpolation of _inData)
                   this case has yet to be implemented
               =3 -perform azimuthal integration or averaging of the 2D intensity distribution _inData and store the resulting 1D distribution in _data and _mesh
                   in that case, the meaning of the subsequent parameters stored in _inPrec is:
                   _inPrec[1] defines whether integration (=0) or averaging (=1) should take place
                   _inPrec[2] defines 1D integration method to be used: =1 means simple integration driven by numbers of points stored in _inPrec[3], =2 means integration driven by relative accuracy specified in _inPrec[3]
                   _inPrec[3] defines number of points (if _inPrec[2] = 1) of relative accuracy (if _inPrec[2] = 2) to be used at the integration
                   _inPrec[4] defines order of interpolation to be used at calculating the values of the distribution _inData out of mesh points
                   _inPrec[5] minimal azimuthal angle for the integration [rad]
                   _inPrec[6] maximal azimuthal angle for the integration [rad]; if _inPrec[6] == _inPrec[5], the integration will be done over 2*Pi 
                   _inPrec[7] horizontal coordinate of center point around which the azimuthal integration should be done
                   _inPrec[8] vertical coordinate of center point around which the azimuthal integration should be done
                   _inPrec[9] list of rectangular areas to be omitted from averaging / integration: [[x0,x_width,y0,y_width],...]
               =4 -fills in ~half of Hermitian Mutual Intensity matrix / data (assuming "normal" data alignment in the complex Hermitian "matrix" E(x,y)*E*(x',y') and filling-out half of it above the diagonal, using complex conjugation)
"""
helpUtiUndFromMagFldTab = """UtiUndFromMagFldTab(_undMagFldC, _inMagFldC, _inPrec)
function attempts to create periodic undulator structure from tabulated magnetic field
:param _undMagFldC: input / output magnetic field container structure (instance of SRWLMagFldC) with undulator structure
       to be set up being allocated in _undMagFldC.arMagFld[0]
:param _inMagFldC: input magnetic field container structure with tabulated field structure to be analyzed;
       the tabulated field structure (instance of SRWLMagFld3D) has to be defined in _inMagFldC.arMagFld[0]
:param _inPrec: input list of precision parameters:
       _inPrec[0]: relative accuracy threshold (nominal value is 0.01)
       _inPrec[1]: maximal number of magnetic field harmonics to attempt to create
       _inPrec[2]: maximal magnetic period length to consider
"""
