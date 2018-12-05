#############################################################################
# SRWLib for Python v 0.14
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import srwlpy as srwl
from array import *
from math import *
from copy import *

import datetime
import json
import random
import sys
import os
import traceback
import uti_math
import errno
import tempfile
import shutil
import time

from srwl_uti_cryst import * 
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
_Light_eV_mu = 1.23984186 #Wavelength <-> Photon Energy conversion constant ([um] <-> [eV])
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
        arHarmLoc = []
        for i in range(_nHarm): arHarm.append(SRWLMagFldH())

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

#****************************************************************************
class SRWLStokes(object):
    """Radiation Stokes Parameters"""
    
    #def __init__(self, _arS0=None, _arS1=None, _arS2=None, _arS3=None, _typeStokes='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0):
    def __init__(self, _arS=None, _typeStokes='f', _eStart=0, _eFin=0, _ne=0, _xStart=0, _xFin=0, _nx=0, _yStart=0, _yFin=0, _ny=0, _mutual=0):
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
        :param _mutual: mutual Stokes components (4*(_ne*_nx*_ny_)^2 values)
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
        if((_arS == 1) and (nProd > 0)):
            self.allocate(_ne, _nx, _ny, _typeStokes, _mutual)          
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
    def allocate(self, _ne, _nx, _ny, _typeStokes='f', _mutual=0):
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
            nTot *= (nTot*2) #OC04052018 (since <E(r1)E(r2)*> may be a complex entity)
        nTot *= 4 #array length of all Stokes components
        #eventually allow for storage of less than 4 Stokes components!

        #DEBUG
        #print('_ne=', _ne, '_nx=', _nx, '_ny=', _ny, 'nTot=', nTot)
        #END DEBUG
        
        self.arS = array(_typeStokes, [0]*nTot)
        self.numTypeStokes = _typeStokes
        self.mesh.ne = _ne
        self.mesh.nx = _nx
        self.mesh.ny = _ny
        self.mutual = _mutual

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

    def avg_update_same_mesh(self, _more_stokes, _iter, _n_stokes_comp=4, _mult=1.):
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
        if(_mult == 1.):
            for ir in range(nStPt):
                self.arS[ir] = (self.arS[ir]*_iter + _more_stokes.arS[ir])/(_iter + 1)
        else:
            for ir in range(nStPt):
                self.arS[ir] = (self.arS[ir]*_iter + _mult*_more_stokes.arS[ir])/(_iter + 1)
            
    def avg_update_interp(self, _more_stokes, _iter, _ord, _n_stokes_comp=4, _mult=1.):
        """ Update this Stokes data structure with new data, contained in the _more_stokes structure, calculated on a different 2D mesh, so that it would represent estimation of average of (_iter + 1) structures
        :param _more_stokes: Stokes data structure to "add" to the estimation of average
        :param _iter: number of Stokes structures already "added" previously
        :param _ord: order of 2D interpolation to use (1- bilinear, ..., 3- bi-cubic)
        :param _n_stokes_comp: number of Stokes components to treat (1 to 4)
        :param _mult: optional multiplier of the _more_stokes
        """

        #DEBUG
        #print('avg_update_interp: iter=', _iter, _mult)
        
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
        xStepWfr = 0
        if(xNpWfr  > 1):
            xStepWfr = (_more_stokes.mesh.xFin - xStartWfr)/(xNpWfr - 1)
        yStartWfr = _more_stokes.mesh.yStart
        yNpWfr = _more_stokes.mesh.ny
        yStepWfr = 0
        if(yNpWfr  > 1):
            yStepWfr = (_more_stokes.mesh.yFin - yStartWfr)/(yNpWfr - 1)
        #DEBUG
        #print('avg_update_interp: iter=', _iter)
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
                        self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)
                        
                        ir += 1
            iOfstSt += nRadWfr        

    def avg_update_interp_mutual(self, _more_stokes, _iter, _n_stokes_comp=4, _mult=1.):
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

                                            #self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)
                                            self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/iter_p_1
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

        orig_iec = int(round(0.5*origNe))
        orig_ixc = int(round(0.5*origNx))
        orig_iyc = int(round(0.5*origNy))
        ofstIc = orig_iec*origPerEp + orig_iec*origPerE + orig_ixc*origPerXp + orig_ixc*origPerX + orig_iyc*origPerYp + orig_iyc*origPerY
        absZerTolI = abs(self.arS[ofstIc])*_rel_zer_tol

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
                                    resDegCohNonRot[resOfst] = absMI/(sqrt(reI1*reI2) + absZerTolI)
                                    resOfst += 1
                        else:
                            ofstMI = ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                            ofstI1 = ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                            ofstI2 = ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY

                            reMI = self.arS[ofstMI]; imMI = self.arS[ofstMI + 1]
                            absMI = sqrt(reMI*reMI + imMI*imMI)
                            reI1 = self.arS[ofstI1]#; imI1 = self.arS[ofstI1 + 1]
                            reI2 = self.arS[ofstI2]#; imI2 = self.arS[ofstI2 + 1]
                            resDegCohNonRot[resOfst] = absMI/(sqrt(reI1*reI2) + absZerTolI)
                            resOfst += 1

        if(not _rot): return resDegCohNonRot

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
        resEstart = origEstart
        resEfin = origEfin
        resXstart = origXstart
        resXfin = origXfin
        resYstart = origYstart
        resYfin = origYfin

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

        resEstep = 0
        if(resNe > 1): resEstep = (resEfin - resEstart)/(resNe - 1)
        resXstep = 0
        if(resNx > 1): resXstep = (resXfin - resXstart)/(resNx - 1)
        resYstep = 0
        if(resNy > 1): resYstep = (resYfin - resYstart)/(resNy - 1)

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
        resYp = resYstart #resY = resYstart
        for iyp in range(resNy): #for iy in range(resNy):
            resY = resYstart #resYp = resYstart
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

                resXp = resXstart #resX = resXstart
                for ixp in range(resNx): #for ix in range(resNx):
                    resX = resXstart #resXp = resXstart
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

                        resEp = resEstart #resE = resEstart
                        for iep in range(resNe): #for ie in range(resNe):
                            resE = resEstart #resEp = resEstart
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
                                else:
                                    resDegCoh[resOfst] = 0

                                resOfst += 1
                                resE += resEstep #resEp += resEstep
                            resEp += resEstep #resE += resEstep
                        resX += resXstep #resXp += resXstep
                    resXp += resXstep #resX += resXstep
                resY += resYstep #resYp += resYstep
            resYp += resYstep #resY += resYstep

        return resDegCoh

    def to_deg_coh_slow(self, _rel_zer_tol=1.e-04, _rot=True):
        """Calculates / "extracts" Degree of Coherence from the Mutual Intensity (first Stokes component, s0)
        :param _rel_zer_tol: relative zero tolerance to use at normalizing (dividing) by the intensity
        :param _rot: rotate or not the degree of coherence data
        :param _n_stokes_comp: number of Stokes components to treat (1 to 4)
        :return: 1D array with (C-aligned) resulting degree of coherence data
        """

        if(self.mutual <= 0): raise Exception("Calculation of Degree of Coherence can not be done from regular Intensity; Mutual Intensity is required.")

        origNe = self.mesh.ne
        origNx = self.mesh.nx
        origNy = self.mesh.ny
        origEstart = self.mesh.eStart
        origEfin = self.mesh.eFin
        origXstart = self.mesh.xStart
        origXfin = self.mesh.xFin
        origYstart = self.mesh.yStart
        origYfin = self.mesh.xFin

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

        origPerEp = 2
        origPerE = origPerEp*origNe
        origPerXp = origPerE*origNe
        origPerX = origPerXp*origNx
        origPerYp = origPerX*origNx
        origPerY = origPerYp*origNy

        orig_iec = int(round(0.5*origNe))
        orig_ixc = int(round(0.5*origNx))
        orig_iyc = int(round(0.5*origNy))
        ofstIc = orig_iec*origPerEp + orig_iec*origPerE + orig_ixc*origPerXp + orig_ixc*origPerX + orig_iyc*origPerYp + orig_iyc*origPerY
        absZerTolI = abs(self.arS[ofstIc])*_rel_zer_tol

        resNe = origNe
        resNx = origNx
        resNy = origNy
        resEstart = origEstart
        resEfin = origEfin
        resXstart = origXstart
        resXfin = origXfin
        resYstart = origYstart
        resYfin = origYfin

        if(_rot):
            sqrt2 = sqrt(2)

            if(resNe > 1): 
                resNe = 2*resNe - 1
                ec = 0.5*(resEstart + resEfin)
                resEstart = ec - sqrt2*(ec - resEstart)
                resEfin = ec + sqrt2*(resEfin - ec)

            if(resNx > 1): 
                resNx = 2*resNx - 1
                xc = 0.5*(resXstart + resXfin)
                resXstart = xc - sqrt2*(xc - resXstart)
                resXfin = xc + sqrt2*(resXfin - xc)

            if(resNy > 1): 
                resNy = 2*resNy - 1
                yc = 0.5*(resYstart + resYfin)
                resYstart = yc - sqrt2*(yc - resYstart)
                resYfin = yc + sqrt2*(resYfin - yc)

        resNp = resNe*resNx*resNy
        resNp *= resNp
        resDegCoh = array('f', [0]*resNp)

        resEstep = 0
        if(resNe > 1): resEstep = (resEfin - resEstart)/(resNe - 1)
        resXstep = 0
        if(resNx > 1): resXstep = (resXfin - resXstart)/(resNx - 1)
        resYstep = 0
        if(resNy > 1): resYstep = (resYfin - resYstart)/(resNy - 1)

        resX1 = 0; resX2 = 0 #just declaring these vars
        resE1 = 0; resE2 = 0 #just declaring these vars

        ie2_0_origPerEp_p_ie1_0_origPerE = 0
        ie2_1_origPerEp_p_ie1_0_origPerE = 0
        ie2_0_origPerEp_p_ie1_1_origPerE = 0
        ie2_1_origPerEp_p_ie1_1_origPerE = 0

        ie1_0_origPerEp_p_ie1_0_origPerE = 0
        ie1_1_origPerEp_p_ie1_0_origPerE = 0
        ie1_0_origPerEp_p_ie1_1_origPerE = 0
        ie1_1_origPerEp_p_ie1_1_origPerE = 0

        ie2_0_origPerEp_p_ie2_0_origPerE = 0
        ie2_1_origPerEp_p_ie2_0_origPerE = 0
        ie2_0_origPerEp_p_ie2_1_origPerE = 0
        ie2_1_origPerEp_p_ie2_1_origPerE = 0

        ix2_0_origPerXp_p_ix1_0_origPerX = 0
        ix2_1_origPerXp_p_ix1_0_origPerX = 0
        ix2_0_origPerXp_p_ix1_1_origPerX = 0
        ix2_1_origPerXp_p_ix1_1_origPerX = 0

        ix1_0_origPerXp_p_ix1_0_origPerX = 0
        ix1_1_origPerXp_p_ix1_0_origPerX = 0
        ix1_0_origPerXp_p_ix1_1_origPerX = 0
        ix1_1_origPerXp_p_ix1_1_origPerX = 0

        ix2_0_origPerXp_p_ix2_0_origPerX = 0
        ix2_1_origPerXp_p_ix2_0_origPerX = 0
        ix2_0_origPerXp_p_ix2_1_origPerX = 0
        ix2_1_origPerXp_p_ix2_1_origPerX = 0

        iy2_0_origPerYp_p_iy1_0_origPerY = 0
        iy2_1_origPerYp_p_iy1_0_origPerY = 0
        iy2_0_origPerYp_p_iy1_1_origPerY = 0
        iy2_1_origPerYp_p_iy1_1_origPerY = 0

        iy1_0_origPerYp_p_iy1_0_origPerY = 0
        iy1_1_origPerYp_p_iy1_0_origPerY = 0
        iy1_0_origPerYp_p_iy1_1_origPerY = 0
        iy1_1_origPerYp_p_iy1_1_origPerY = 0

        iy2_0_origPerYp_p_iy2_0_origPerY = 0
        iy2_1_origPerYp_p_iy2_0_origPerY = 0
        iy2_0_origPerYp_p_iy2_1_origPerY = 0
        iy2_1_origPerYp_p_iy2_1_origPerY = 0

        resOfst = 0
        resY = resYstart
        for iy in range(resNy):
            resYp = resYstart
            for iyp in range(resNy):
                isInRangeY = True
                if(_rot and (resNy > 1)):
                    resY1 = resY + resYp
                    resY2 = resY - resYp
                    #resY1 = resY - resYp
                    #resY2 = resY + resYp
                    if((resY1 < origYstart) or (resY1 > origYfin) or (resY2 < origYstart) or (resY2 > origYfin)): isInRangeY = False

                iy1_0 = 0; iy1_1 = 0; ry1 = 0
                iy2_0 = 0; iy2_1 = 0; ry2 = 0
                if(_rot):
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

                        iy2_0_origPerYp = iy2_0*origPerYp
                        iy2_1_origPerYp = iy2_1*origPerYp
                        iy1_0_origPerY = iy1_0*origPerY
                        iy1_1_origPerY = iy1_1*origPerY
                        iy2_0_origPerYp_p_iy1_0_origPerY = iy2_0_origPerYp + iy1_0_origPerY
                        iy2_1_origPerYp_p_iy1_0_origPerY = iy2_1_origPerYp + iy1_0_origPerY
                        iy2_0_origPerYp_p_iy1_1_origPerY = iy2_0_origPerYp + iy1_1_origPerY
                        iy2_1_origPerYp_p_iy1_1_origPerY = iy2_1_origPerYp + iy1_1_origPerY
                        iy1_0_origPerYp = iy1_0*origPerYp
                        iy1_1_origPerYp = iy1_1*origPerYp
                        iy1_0_origPerYp_p_iy1_0_origPerY = iy1_0_origPerYp + iy1_0_origPerY
                        iy1_1_origPerYp_p_iy1_0_origPerY = iy1_1_origPerYp + iy1_0_origPerY
                        iy1_0_origPerYp_p_iy1_1_origPerY = iy1_0_origPerYp + iy1_1_origPerY
                        iy1_1_origPerYp_p_iy1_1_origPerY = iy1_1_origPerYp + iy1_1_origPerY
                        iy2_0_origPerY = iy2_0*origPerY
                        iy2_1_origPerY = iy2_1*origPerY
                        iy2_0_origPerYp_p_iy2_0_origPerY = iy2_0_origPerYp + iy2_0_origPerY
                        iy2_1_origPerYp_p_iy2_0_origPerY = iy2_1_origPerYp + iy2_0_origPerY
                        iy2_0_origPerYp_p_iy2_1_origPerY = iy2_0_origPerYp + iy2_1_origPerY
                        iy2_1_origPerYp_p_iy2_1_origPerY = iy2_1_origPerYp + iy2_1_origPerY
                else:
                    iy1_0 = iy
                    iy2_0 = iyp
                    iy2_0_origPerYp = iy2_0*origPerYp
                    iy1_0_origPerY = iy1_0*origPerY
                    iy2_0_origPerYp_p_iy1_0_origPerY = iy2_0_origPerYp + iy1_0_origPerY
                    iy1_0_origPerYp_p_iy1_0_origPerY = iy1_0*origPerYp + iy1_0_origPerY
                    iy2_0_origPerYp_p_iy2_0_origPerY = iy2_0_origPerYp + iy2_0*origPerY

                resX = resXstart
                for ix in range(resNx):
                    resXp = resXstart
                    for ixp in range(resNx):
                        isInRangeX = True
                        if(isInRangeY):
                            if(_rot and (resNx > 1)):
                                resX1 = resX + resXp
                                resX2 = resX - resXp
                                #resX1 = resX - resXp
                                #resX2 = resX + resXp
                                if((resX1 < origXstart) or (resX1 > origXfin) or (resX2 < origXstart) or (resX2 > origXfin)): isInRangeX = False
                        else:
                            isInRangeX = False

                        ix1_0 = 0; ix1_1 = 0; rx1 = 0
                        ix2_0 = 0; ix2_1 = 0; rx2 = 0
                        if(_rot):
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

                                ix2_0_origPerXp = ix2_0*origPerXp
                                ix2_1_origPerXp = ix2_1*origPerXp
                                ix1_0_origPerX = ix1_0*origPerX
                                ix1_1_origPerX = ix1_1*origPerX
                                ix2_0_origPerXp_p_ix1_0_origPerX = ix2_0_origPerXp + ix1_0_origPerX
                                ix2_1_origPerXp_p_ix1_0_origPerX = ix2_1_origPerXp + ix1_0_origPerX
                                ix2_0_origPerXp_p_ix1_1_origPerX = ix2_0_origPerXp + ix1_1_origPerX
                                ix2_1_origPerXp_p_ix1_1_origPerX = ix2_1_origPerXp + ix1_1_origPerX
                                ix1_0_origPerXp = ix1_0*origPerXp
                                ix1_1_origPerXp = ix1_1*origPerXp
                                ix1_0_origPerXp_p_ix1_0_origPerX = ix1_0_origPerXp + ix1_0_origPerX
                                ix1_1_origPerXp_p_ix1_0_origPerX = ix1_1_origPerXp + ix1_0_origPerX
                                ix1_0_origPerXp_p_ix1_1_origPerX = ix1_0_origPerXp + ix1_1_origPerX
                                ix1_1_origPerXp_p_ix1_1_origPerX = ix1_1_origPerXp + ix1_1_origPerX
                                ix2_0_origPerX = ix2_0*origPerX
                                ix2_1_origPerX = ix2_1*origPerX
                                ix2_0_origPerXp_p_ix2_0_origPerX = ix2_0_origPerXp + ix2_0_origPerX
                                ix2_1_origPerXp_p_ix2_0_origPerX = ix2_1_origPerXp + ix2_0_origPerX
                                ix2_0_origPerXp_p_ix2_1_origPerX = ix2_0_origPerXp + ix2_1_origPerX
                                ix2_1_origPerXp_p_ix2_1_origPerX = ix2_1_origPerXp + ix2_1_origPerX
                        else:
                            ix1_0 = ix
                            ix2_0 = ixp
                            ix2_0_origPerXp = ix2_0*origPerXp
                            ix1_0_origPerX = ix1_0*origPerX
                            ix2_0_origPerXp_p_ix1_0_origPerX = ix2_0_origPerXp + ix1_0_origPerX
                            ix1_0_origPerXp_p_ix1_0_origPerX = ix1_0*origPerXp + ix1_0_origPerX
                            ix2_0_origPerXp_p_ix2_0_origPerX = ix2_0_origPerXp + ix2_0*origPerX

                        resE = resEstart
                        for ie in range(resNe):
                            resEp = resEstart
                            for iep in range(resNe):
                                isInRangeE = True
                                if(isInRangeX):
                                    if(_rot and (resNe > 1)):
                                        resE1 = resE + resEp
                                        resE2 = resE - resEp
                                        #resE1 = resE - resEp
                                        #resE2 = resE + resEp
                                        if((resE1 < origEstart) or (resE1 > origEfin) or (resE2 < origEstart) or (resE2 > origEfin)): isInRangeE = False
                                else:
                                    isInRangeE = False

                                ie1_0 = 0; ie1_1 = 0; re1 = 0
                                ie2_0 = 0; ie2_1 = 0; re2 = 0
                                if(_rot):
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

                                        ie2_0_origPerEp = ie2_0*origPerEp
                                        ie2_1_origPerEp = ie2_1*origPerEp
                                        ie1_0_origPerE = ie1_0*origPerE
                                        ie1_1_origPerE = ie1_1*origPerE
                                        ie2_0_origPerEp_p_ie1_0_origPerE = ie2_0_origPerEp + ie1_0_origPerE
                                        ie2_1_origPerEp_p_ie1_0_origPerE = ie2_1_origPerEp + ie1_0_origPerE
                                        ie2_0_origPerEp_p_ie1_1_origPerE = ie2_0_origPerEp + ie1_1_origPerE
                                        ie2_1_origPerEp_p_ie1_1_origPerE = ie2_1_origPerEp + ie1_1_origPerE
                                        ie1_0_origPerEp = ie1_0*origPerEp
                                        ie1_1_origPerEp = ie1_1*origPerEp
                                        ie1_0_origPerEp_p_ie1_0_origPerE = ie1_0_origPerEp + ie1_0_origPerE
                                        ie1_1_origPerEp_p_ie1_0_origPerE = ie1_1_origPerEp + ie1_0_origPerE
                                        ie1_0_origPerEp_p_ie1_1_origPerE = ie1_0_origPerEp + ie1_1_origPerE
                                        ie1_1_origPerEp_p_ie1_1_origPerE = ie1_1_origPerEp + ie1_1_origPerE
                                        ie2_0_origPerE = ie2_0*origPerE
                                        ie2_1_origPerE = ie2_1*origPerE
                                        ie2_0_origPerEp_p_ie2_0_origPerE = ie2_0_origPerEp + ie2_0_origPerE
                                        ie2_1_origPerEp_p_ie2_0_origPerE = ie2_1_origPerEp + ie2_0_origPerE
                                        ie2_0_origPerEp_p_ie2_1_origPerE = ie2_0_origPerEp + ie2_1_origPerE
                                        ie2_1_origPerEp_p_ie2_1_origPerE = ie2_1_origPerEp + ie2_1_origPerE

                                    if(isInRangeE):
                                        #--------------MI
                                        ofst000000 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst100000 = ie2_1_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst010000 = ie2_0_origPerEp_p_ie1_1_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst001000 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_1_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst000100 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_1_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst000010 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_1_origPerYp_p_iy1_0_origPerY
                                        ofst000001 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_1_origPerY

                                        ofst110000 = ie2_1_origPerEp_p_ie1_1_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst101000 = ie2_1_origPerEp_p_ie1_0_origPerE + ix2_1_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst100100 = ie2_1_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_1_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst100010 = ie2_1_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_1_origPerYp_p_iy1_0_origPerY
                                        ofst100001 = ie2_1_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_1_origPerY
                                    
                                        ofst011000 = ie2_0_origPerEp_p_ie1_1_origPerE + ix2_1_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst010100 = ie2_0_origPerEp_p_ie1_1_origPerE + ix2_0_origPerXp_p_ix1_1_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst010010 = ie2_0_origPerEp_p_ie1_1_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_1_origPerYp_p_iy1_0_origPerY
                                        ofst010001 = ie2_0_origPerEp_p_ie1_1_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_1_origPerY
                                   
                                        ofst001100 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_1_origPerXp_p_ix1_1_origPerX + iy2_0_origPerYp_p_iy1_0_origPerY
                                        ofst001010 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_1_origPerXp_p_ix1_0_origPerX + iy2_1_origPerYp_p_iy1_0_origPerY
                                        ofst001001 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_1_origPerXp_p_ix1_0_origPerX + iy2_0_origPerYp_p_iy1_1_origPerY

                                        ofst000110 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_1_origPerX + iy2_1_origPerYp_p_iy1_0_origPerY
                                        ofst000101 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_1_origPerX + iy2_0_origPerYp_p_iy1_1_origPerY

                                        ofst000011 = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0_origPerXp_p_ix1_0_origPerX + iy2_1_origPerYp_p_iy1_1_origPerY
                                    
                                        re_a000000 = _more_stokes.arS[ofst000000]; ofst000000 += 1; im_a000000 = _more_stokes.arS[ofst000000]
                                        re_f100000 = _more_stokes.arS[ofst100000]; ofst100000 += 1; im_f100000 = _more_stokes.arS[ofst100000]
                                        re_f010000 = _more_stokes.arS[ofst010000]; ofst010000 += 1; im_f010000 = _more_stokes.arS[ofst010000]
                                        re_f001000 = _more_stokes.arS[ofst001000]; ofst001000 += 1; im_f001000 = _more_stokes.arS[ofst001000]
                                        re_f000100 = _more_stokes.arS[ofst000100]; ofst000100 += 1; im_f000100 = _more_stokes.arS[ofst000100]
                                        re_f000010 = _more_stokes.arS[ofst000010]; ofst000010 += 1; im_f000010 = _more_stokes.arS[ofst000010]
                                        re_f000001 = _more_stokes.arS[ofst000001]; ofst000001 += 1; im_f000001 = _more_stokes.arS[ofst000001]

                                        re_f110000 = _more_stokes.arS[ofst110000]; ofst110000 += 1; im_f110000 = _more_stokes.arS[ofst110000]
                                        re_f101000 = _more_stokes.arS[ofst101000]; ofst101000 += 1; im_f101000 = _more_stokes.arS[ofst101000]
                                        re_f100100 = _more_stokes.arS[ofst100100]; ofst100100 += 1; im_f100100 = _more_stokes.arS[ofst100100]
                                        re_f100010 = _more_stokes.arS[ofst100010]; ofst100010 += 1; im_f100010 = _more_stokes.arS[ofst100010]
                                        re_f100001 = _more_stokes.arS[ofst100001]; ofst100001 += 1; im_f100001 = _more_stokes.arS[ofst100001]

                                        re_f011000 = _more_stokes.arS[ofst011000]; ofst011000 += 1; im_f011000 = _more_stokes.arS[ofst011000]
                                        re_f010100 = _more_stokes.arS[ofst010100]; ofst010100 += 1; im_f010100 = _more_stokes.arS[ofst010100]
                                        re_f010010 = _more_stokes.arS[ofst010010]; ofst010010 += 1; im_f010010 = _more_stokes.arS[ofst010010]
                                        re_f010001 = _more_stokes.arS[ofst010001]; ofst010001 += 1; im_f010001 = _more_stokes.arS[ofst010001]

                                        re_f001100 = _more_stokes.arS[ofst001100]; ofst001100 += 1; im_f001100 = _more_stokes.arS[ofst001100]
                                        re_f001010 = _more_stokes.arS[ofst001010]; ofst001010 += 1; im_f001010 = _more_stokes.arS[ofst001010]
                                        re_f001001 = _more_stokes.arS[ofst001001]; ofst001001 += 1; im_f001001 = _more_stokes.arS[ofst001001]

                                        re_f000110 = _more_stokes.arS[ofst000110]; ofst000110 += 1; im_f000110 = _more_stokes.arS[ofst000110]
                                        re_f000101 = _more_stokes.arS[ofst000101]; ofst000101 += 1; im_f000101 = _more_stokes.arS[ofst000101]
                                        
                                        re_f000011 = _more_stokes.arS[ofst000011]; ofst000011 += 1; im_f000011 = _more_stokes.arS[ofst000011]

                                        re_a100000 = re_f100000 - re_a000000; im_a100000 = im_f100000 - im_a000000
                                        re_a010000 = re_f010000 - re_a000000; im_a010000 = im_f010000 - im_a000000
                                        re_a001000 = re_f001000 - re_a000000; im_a001000 = im_f001000 - im_a000000
                                        re_a000100 = re_f000100 - re_a000000; im_a000100 = im_f000100 - im_a000000
                                        re_a000010 = re_f000010 - re_a000000; im_a000010 = im_f000010 - im_a000000
                                        re_a000001 = re_f000001 - re_a000000; im_a000001 = im_f000001 - im_a000000
                                        re_a110000 = re_a000000 - re_f010000 - re_f100000 + re_f110000; im_a110000 = im_a000000 - im_f010000 - im_f100000 + im_f110000
                                        re_a101000 = re_a000000 - re_f001000 - re_f100000 + re_f101000; im_a101000 = im_a000000 - im_f001000 - im_f100000 + im_f101000
                                        re_a100100 = re_a000000 - re_f000100 - re_f100000 + re_f100100; im_a100100 = im_a000000 - im_f000100 - im_f100000 + im_f100100
                                        re_a100010 = re_a000000 - re_f000010 - re_f100000 + re_f100010; im_a100010 = im_a000000 - im_f000010 - im_f100000 + im_f100010
                                        re_a100001 = re_a000000 - re_f000001 - re_f100000 + re_f100001; im_a100001 = im_a000000 - im_f000001 - im_f100000 + im_f100001
                                        re_a011000 = re_a000000 - re_f001000 - re_f010000 + re_f011000; im_a011000 = im_a000000 - im_f001000 - im_f010000 + im_f011000
                                        re_a010100 = re_a000000 - re_f000100 - re_f010000 + re_f010100; im_a010100 = im_a000000 - im_f000100 - im_f010000 + im_f010100
                                        re_a010010 = re_a000000 - re_f000010 - re_f010000 + re_f010010; im_a010010 = im_a000000 - im_f000010 - im_f010000 + im_f010010
                                        re_a010001 = re_a000000 - re_f000001 - re_f010000 + re_f010001; im_a010001 = im_a000000 - im_f000001 - im_f010000 + im_f010001
                                        re_a001100 = re_a000000 - re_f000100 - re_f001000 + re_f001100; im_a001100 = im_a000000 - im_f000100 - im_f001000 + im_f001100
                                        re_a001010 = re_a000000 - re_f000010 - re_f001000 + re_f001010; im_a001010 = im_a000000 - im_f000010 - im_f001000 + im_f001010
                                        re_a001001 = re_a000000 - re_f000001 - re_f001000 + re_f001001; im_a001001 = im_a000000 - im_f000001 - im_f001000 + im_f001001
                                        re_a000110 = re_a000000 - re_f000010 - re_f000100 + re_f000110; im_a000110 = im_a000000 - im_f000010 - im_f000100 + im_f000110
                                        re_a000101 = re_a000000 - re_f000001 - re_f000100 + re_f000101; im_a000101 = im_a000000 - im_f000001 - im_f000100 + im_f000101
                                        re_a000011 = re_a000000 - re_f000001 - re_f000010 + re_f000011; im_a000011 = im_a000000 - im_f000001 - im_f000010 + im_f000011

                                        reMI = (re_a100000 + re_a110000*re1 + re_a101000*rx2 + re_a100100*rx1 + re_a100010*ry2 + re_a100001*ry1)*re2
                                        reMI += (re_a010000 + re_a011000*rx2 + re_a010100*rx1 + re_a010010*ry2 + re_a010001*ry1)*re1
                                        reMI += (re_a001000 + re_a001100*rx1 + re_a001010*ry2 + re_a001001*ry1)*rx2
                                        reMI += (re_a000100 + re_a000110*ry2 + re_a000101*ry1)*rx1 + (re_a000010 + re_a000011*ry1)*ry2 + re_a000001*ry1 + re_a000000
                                        imMI = (im_a100000 + im_a110000*re1 + im_a101000*rx2 + im_a100100*rx1 + im_a100010*ry2 + im_a100001*ry1)*re2
                                        imMI += (im_a010000 + im_a011000*rx2 + im_a010100*rx1 + im_a010010*ry2 + im_a010001*ry1)*re1
                                        imMI += (im_a001000 + im_a001100*rx1 + im_a001010*ry2 + im_a001001*ry1)*rx2
                                        imMI += (im_a000100 + im_a000110*ry2 + im_a000101*ry1)*rx1 + (im_a000010 + im_a000011*ry1)*ry2 + im_a000001*ry1 + im_a000000
                                        absMI = sqrt(reMI*reMI + imMI*imMI)

                                        #--------------I1
                                        ofst000000 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst100000 = ie1_1_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst010000 = ie1_0_origPerEp_p_ie1_1_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst001000 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_1_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst000100 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_1_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst000010 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_1_origPerYp_p_iy1_0_origPerY
                                        ofst000001 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_1_origPerY

                                        ofst110000 = ie1_1_origPerEp_p_ie1_1_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst101000 = ie1_1_origPerEp_p_ie1_0_origPerE + ix1_1_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst100100 = ie1_1_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_1_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst100010 = ie1_1_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_1_origPerYp_p_iy1_0_origPerY
                                        ofst100001 = ie1_1_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_1_origPerY
                                    
                                        ofst011000 = ie1_0_origPerEp_p_ie1_1_origPerE + ix1_1_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst010100 = ie1_0_origPerEp_p_ie1_1_origPerE + ix1_0_origPerXp_p_ix1_1_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst010010 = ie1_0_origPerEp_p_ie1_1_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_1_origPerYp_p_iy1_0_origPerY
                                        ofst010001 = ie1_0_origPerEp_p_ie1_1_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_1_origPerY
                                   
                                        ofst001100 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_1_origPerXp_p_ix1_1_origPerX + iy1_0_origPerYp_p_iy1_0_origPerY
                                        ofst001010 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_1_origPerXp_p_ix1_0_origPerX + iy1_1_origPerYp_p_iy1_0_origPerY
                                        ofst001001 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_1_origPerXp_p_ix1_0_origPerX + iy1_0_origPerYp_p_iy1_1_origPerY

                                        ofst000110 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_1_origPerX + iy1_1_origPerYp_p_iy1_0_origPerY
                                        ofst000101 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_1_origPerX + iy1_0_origPerYp_p_iy1_1_origPerY

                                        ofst000011 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0_origPerXp_p_ix1_0_origPerX + iy1_1_origPerYp_p_iy1_1_origPerY
                                    
                                        re_a000000 = _more_stokes.arS[ofst000000]#; ofst000000 += 1; im_a000000 = _more_stokes.arS[ofst000000]
                                        re_f100000 = _more_stokes.arS[ofst100000]#; ofst100000 += 1; im_f100000 = _more_stokes.arS[ofst100000]
                                        re_f010000 = _more_stokes.arS[ofst010000]#; ofst010000 += 1; im_f010000 = _more_stokes.arS[ofst010000]
                                        re_f001000 = _more_stokes.arS[ofst001000]#; ofst001000 += 1; im_f001000 = _more_stokes.arS[ofst001000]
                                        re_f000100 = _more_stokes.arS[ofst000100]#; ofst000100 += 1; im_f000100 = _more_stokes.arS[ofst000100]
                                        re_f000010 = _more_stokes.arS[ofst000010]#; ofst000010 += 1; im_f000010 = _more_stokes.arS[ofst000010]
                                        re_f000001 = _more_stokes.arS[ofst000001]#; ofst000001 += 1; im_f000001 = _more_stokes.arS[ofst000001]

                                        re_f110000 = _more_stokes.arS[ofst110000]#; ofst110000 += 1; im_f110000 = _more_stokes.arS[ofst110000]
                                        re_f101000 = _more_stokes.arS[ofst101000]#; ofst101000 += 1; im_f101000 = _more_stokes.arS[ofst101000]
                                        re_f100100 = _more_stokes.arS[ofst100100]#; ofst100100 += 1; im_f100100 = _more_stokes.arS[ofst100100]
                                        re_f100010 = _more_stokes.arS[ofst100010]#; ofst100010 += 1; im_f100010 = _more_stokes.arS[ofst100010]
                                        re_f100001 = _more_stokes.arS[ofst100001]#; ofst100001 += 1; im_f100001 = _more_stokes.arS[ofst100001]

                                        re_f011000 = _more_stokes.arS[ofst011000]#; ofst011000 += 1; im_f011000 = _more_stokes.arS[ofst011000]
                                        re_f010100 = _more_stokes.arS[ofst010100]#; ofst010100 += 1; im_f010100 = _more_stokes.arS[ofst010100]
                                        re_f010010 = _more_stokes.arS[ofst010010]#; ofst010010 += 1; im_f010010 = _more_stokes.arS[ofst010010]
                                        re_f010001 = _more_stokes.arS[ofst010001]#; ofst010001 += 1; im_f010001 = _more_stokes.arS[ofst010001]

                                        re_f001100 = _more_stokes.arS[ofst001100]#; ofst001100 += 1; im_f001100 = _more_stokes.arS[ofst001100]
                                        re_f001010 = _more_stokes.arS[ofst001010]#; ofst001010 += 1; im_f001010 = _more_stokes.arS[ofst001010]
                                        re_f001001 = _more_stokes.arS[ofst001001]#; ofst001001 += 1; im_f001001 = _more_stokes.arS[ofst001001]

                                        re_f000110 = _more_stokes.arS[ofst000110]#; ofst000110 += 1; im_f000110 = _more_stokes.arS[ofst000110]
                                        re_f000101 = _more_stokes.arS[ofst000101]#; ofst000101 += 1; im_f000101 = _more_stokes.arS[ofst000101]
                                        
                                        re_f000011 = _more_stokes.arS[ofst000011]#; ofst000011 += 1; im_f000011 = _more_stokes.arS[ofst000011]

                                        re_a100000 = re_f100000 - re_a000000#; im_a100000 = im_f100000 - im_a000000
                                        re_a010000 = re_f010000 - re_a000000#; im_a010000 = im_f010000 - im_a000000
                                        re_a001000 = re_f001000 - re_a000000#; im_a001000 = im_f001000 - im_a000000
                                        re_a000100 = re_f000100 - re_a000000#; im_a000100 = im_f000100 - im_a000000
                                        re_a000010 = re_f000010 - re_a000000#; im_a000010 = im_f000010 - im_a000000
                                        re_a000001 = re_f000001 - re_a000000#; im_a000001 = im_f000001 - im_a000000
                                        re_a110000 = re_a000000 - re_f010000 - re_f100000 + re_f110000#; im_a110000 = im_a000000 - im_f010000 - im_f100000 + im_f110000
                                        re_a101000 = re_a000000 - re_f001000 - re_f100000 + re_f101000#; im_a101000 = im_a000000 - im_f001000 - im_f100000 + im_f101000
                                        re_a100100 = re_a000000 - re_f000100 - re_f100000 + re_f100100#; im_a100100 = im_a000000 - im_f000100 - im_f100000 + im_f100100
                                        re_a100010 = re_a000000 - re_f000010 - re_f100000 + re_f100010#; im_a100010 = im_a000000 - im_f000010 - im_f100000 + im_f100010
                                        re_a100001 = re_a000000 - re_f000001 - re_f100000 + re_f100001#; im_a100001 = im_a000000 - im_f000001 - im_f100000 + im_f100001
                                        re_a011000 = re_a000000 - re_f001000 - re_f010000 + re_f011000#; im_a011000 = im_a000000 - im_f001000 - im_f010000 + im_f011000
                                        re_a010100 = re_a000000 - re_f000100 - re_f010000 + re_f010100#; im_a010100 = im_a000000 - im_f000100 - im_f010000 + im_f010100
                                        re_a010010 = re_a000000 - re_f000010 - re_f010000 + re_f010010#; im_a010010 = im_a000000 - im_f000010 - im_f010000 + im_f010010
                                        re_a010001 = re_a000000 - re_f000001 - re_f010000 + re_f010001#; im_a010001 = im_a000000 - im_f000001 - im_f010000 + im_f010001
                                        re_a001100 = re_a000000 - re_f000100 - re_f001000 + re_f001100#; im_a001100 = im_a000000 - im_f000100 - im_f001000 + im_f001100
                                        re_a001010 = re_a000000 - re_f000010 - re_f001000 + re_f001010#; im_a001010 = im_a000000 - im_f000010 - im_f001000 + im_f001010
                                        re_a001001 = re_a000000 - re_f000001 - re_f001000 + re_f001001#; im_a001001 = im_a000000 - im_f000001 - im_f001000 + im_f001001
                                        re_a000110 = re_a000000 - re_f000010 - re_f000100 + re_f000110#; im_a000110 = im_a000000 - im_f000010 - im_f000100 + im_f000110
                                        re_a000101 = re_a000000 - re_f000001 - re_f000100 + re_f000101#; im_a000101 = im_a000000 - im_f000001 - im_f000100 + im_f000101
                                        re_a000011 = re_a000000 - re_f000001 - re_f000010 + re_f000011#; im_a000011 = im_a000000 - im_f000001 - im_f000010 + im_f000011

                                        reI1 = (re_a100000 + re_a110000*re1 + (re_a101000 + re_a100100)*rx1 + (re_a100010 + re_a100001)*ry1)*re1
                                        reI1 += (re_a010000 + (re_a011000 + re_a010100)*rx1 + (re_a010010 + re_a010001)*ry1)*re1
                                        reI1 += (re_a001000 + re_a001100*rx1 + (re_a001010 + re_a001001)*ry1)*rx1
                                        reI1 += (re_a000100 + (re_a000110 + re_a000101)*ry1)*rx1 + (re_a000010 + re_a000011*ry1)*ry1 + re_a000001*ry1 + re_a000000
                                        #imI1 = (im_a100000 + im_a110000*re1 + (im_a101000 + im_a100100)*rx1 + (im_a100010 + im_a100001)*ry1)*re1
                                        #imI1 += (im_a010000 + (im_a011000 + im_a010100)*rx1 + (im_a010010 + im_a010001)*ry1)*re1
                                        #imI1 += (im_a001000 + im_a001100*rx1 + (im_a001010 + im_a001001)*ry1)*rx1
                                        #imI1 += (im_a000100 + (im_a000110 + im_a000101)*ry1)*rx1 + (im_a000010 + im_a000011*ry1)*ry1 + im_a000001*ry1 + im_a000000

                                        #--------------I2
                                        ofst000000 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst100000 = ie2_1_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst010000 = ie2_0_origPerEp_p_ie2_1_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst001000 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_1_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst000100 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_1_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst000010 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_1_origPerYp_p_iy2_0_origPerY
                                        ofst000001 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_1_origPerY

                                        ofst110000 = ie2_1_origPerEp_p_ie2_1_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst101000 = ie2_1_origPerEp_p_ie2_0_origPerE + ix2_1_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst100100 = ie2_1_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_1_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst100010 = ie2_1_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_1_origPerYp_p_iy2_0_origPerY
                                        ofst100001 = ie2_1_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_1_origPerY
                                    
                                        ofst011000 = ie2_0_origPerEp_p_ie2_1_origPerE + ix2_1_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst010100 = ie2_0_origPerEp_p_ie2_1_origPerE + ix2_0_origPerXp_p_ix2_1_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst010010 = ie2_0_origPerEp_p_ie2_1_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_1_origPerYp_p_iy2_0_origPerY
                                        ofst010001 = ie2_0_origPerEp_p_ie2_1_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_1_origPerY
                                   
                                        ofst001100 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_1_origPerXp_p_ix2_1_origPerX + iy2_0_origPerYp_p_iy2_0_origPerY
                                        ofst001010 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_1_origPerXp_p_ix2_0_origPerX + iy2_1_origPerYp_p_iy2_0_origPerY
                                        ofst001001 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_1_origPerXp_p_ix2_0_origPerX + iy2_0_origPerYp_p_iy2_1_origPerY

                                        ofst000110 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_1_origPerX + iy2_1_origPerYp_p_iy2_0_origPerY
                                        ofst000101 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_1_origPerX + iy2_0_origPerYp_p_iy2_1_origPerY

                                        ofst000011 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0_origPerXp_p_ix2_0_origPerX + iy2_1_origPerYp_p_iy2_1_origPerY
                                    
                                        re_a000000 = _more_stokes.arS[ofst000000]#; ofst000000 += 1; im_a000000 = _more_stokes.arS[ofst000000]
                                        re_f100000 = _more_stokes.arS[ofst100000]#; ofst100000 += 1; im_f100000 = _more_stokes.arS[ofst100000]
                                        re_f010000 = _more_stokes.arS[ofst010000]#; ofst010000 += 1; im_f010000 = _more_stokes.arS[ofst010000]
                                        re_f001000 = _more_stokes.arS[ofst001000]#; ofst001000 += 1; im_f001000 = _more_stokes.arS[ofst001000]
                                        re_f000100 = _more_stokes.arS[ofst000100]#; ofst000100 += 1; im_f000100 = _more_stokes.arS[ofst000100]
                                        re_f000010 = _more_stokes.arS[ofst000010]#; ofst000010 += 1; im_f000010 = _more_stokes.arS[ofst000010]
                                        re_f000001 = _more_stokes.arS[ofst000001]#; ofst000001 += 1; im_f000001 = _more_stokes.arS[ofst000001]

                                        re_f110000 = _more_stokes.arS[ofst110000]#; ofst110000 += 1; im_f110000 = _more_stokes.arS[ofst110000]
                                        re_f101000 = _more_stokes.arS[ofst101000]#; ofst101000 += 1; im_f101000 = _more_stokes.arS[ofst101000]
                                        re_f100100 = _more_stokes.arS[ofst100100]#; ofst100100 += 1; im_f100100 = _more_stokes.arS[ofst100100]
                                        re_f100010 = _more_stokes.arS[ofst100010]#; ofst100010 += 1; im_f100010 = _more_stokes.arS[ofst100010]
                                        re_f100001 = _more_stokes.arS[ofst100001]#; ofst100001 += 1; im_f100001 = _more_stokes.arS[ofst100001]

                                        re_f011000 = _more_stokes.arS[ofst011000]#; ofst011000 += 1; im_f011000 = _more_stokes.arS[ofst011000]
                                        re_f010100 = _more_stokes.arS[ofst010100]#; ofst010100 += 1; im_f010100 = _more_stokes.arS[ofst010100]
                                        re_f010010 = _more_stokes.arS[ofst010010]#; ofst010010 += 1; im_f010010 = _more_stokes.arS[ofst010010]
                                        re_f010001 = _more_stokes.arS[ofst010001]#; ofst010001 += 1; im_f010001 = _more_stokes.arS[ofst010001]

                                        re_f001100 = _more_stokes.arS[ofst001100]#; ofst001100 += 1; im_f001100 = _more_stokes.arS[ofst001100]
                                        re_f001010 = _more_stokes.arS[ofst001010]#; ofst001010 += 1; im_f001010 = _more_stokes.arS[ofst001010]
                                        re_f001001 = _more_stokes.arS[ofst001001]#; ofst001001 += 1; im_f001001 = _more_stokes.arS[ofst001001]

                                        re_f000110 = _more_stokes.arS[ofst000110]#; ofst000110 += 1; im_f000110 = _more_stokes.arS[ofst000110]
                                        re_f000101 = _more_stokes.arS[ofst000101]#; ofst000101 += 1; im_f000101 = _more_stokes.arS[ofst000101]
                                        
                                        re_f000011 = _more_stokes.arS[ofst000011]#; ofst000011 += 1; im_f000011 = _more_stokes.arS[ofst000011]

                                        re_a100000 = re_f100000 - re_a000000#; im_a100000 = im_f100000 - im_a000000
                                        re_a010000 = re_f010000 - re_a000000#; im_a010000 = im_f010000 - im_a000000
                                        re_a001000 = re_f001000 - re_a000000#; im_a001000 = im_f001000 - im_a000000
                                        re_a000100 = re_f000100 - re_a000000#; im_a000100 = im_f000100 - im_a000000
                                        re_a000010 = re_f000010 - re_a000000#; im_a000010 = im_f000010 - im_a000000
                                        re_a000001 = re_f000001 - re_a000000#; im_a000001 = im_f000001 - im_a000000
                                        re_a110000 = re_a000000 - re_f010000 - re_f100000 + re_f110000#; im_a110000 = im_a000000 - im_f010000 - im_f100000 + im_f110000
                                        re_a101000 = re_a000000 - re_f001000 - re_f100000 + re_f101000#; im_a101000 = im_a000000 - im_f001000 - im_f100000 + im_f101000
                                        re_a100100 = re_a000000 - re_f000100 - re_f100000 + re_f100100#; im_a100100 = im_a000000 - im_f000100 - im_f100000 + im_f100100
                                        re_a100010 = re_a000000 - re_f000010 - re_f100000 + re_f100010#; im_a100010 = im_a000000 - im_f000010 - im_f100000 + im_f100010
                                        re_a100001 = re_a000000 - re_f000001 - re_f100000 + re_f100001#; im_a100001 = im_a000000 - im_f000001 - im_f100000 + im_f100001
                                        re_a011000 = re_a000000 - re_f001000 - re_f010000 + re_f011000#; im_a011000 = im_a000000 - im_f001000 - im_f010000 + im_f011000
                                        re_a010100 = re_a000000 - re_f000100 - re_f010000 + re_f010100#; im_a010100 = im_a000000 - im_f000100 - im_f010000 + im_f010100
                                        re_a010010 = re_a000000 - re_f000010 - re_f010000 + re_f010010#; im_a010010 = im_a000000 - im_f000010 - im_f010000 + im_f010010
                                        re_a010001 = re_a000000 - re_f000001 - re_f010000 + re_f010001#; im_a010001 = im_a000000 - im_f000001 - im_f010000 + im_f010001
                                        re_a001100 = re_a000000 - re_f000100 - re_f001000 + re_f001100#; im_a001100 = im_a000000 - im_f000100 - im_f001000 + im_f001100
                                        re_a001010 = re_a000000 - re_f000010 - re_f001000 + re_f001010#; im_a001010 = im_a000000 - im_f000010 - im_f001000 + im_f001010
                                        re_a001001 = re_a000000 - re_f000001 - re_f001000 + re_f001001#; im_a001001 = im_a000000 - im_f000001 - im_f001000 + im_f001001
                                        re_a000110 = re_a000000 - re_f000010 - re_f000100 + re_f000110#; im_a000110 = im_a000000 - im_f000010 - im_f000100 + im_f000110
                                        re_a000101 = re_a000000 - re_f000001 - re_f000100 + re_f000101#; im_a000101 = im_a000000 - im_f000001 - im_f000100 + im_f000101
                                        re_a000011 = re_a000000 - re_f000001 - re_f000010 + re_f000011#; im_a000011 = im_a000000 - im_f000001 - im_f000010 + im_f000011

                                        reI2 = (re_a100000 + re_a110000*re2 + (re_a101000 + re_a100100)*rx2 + (re_a100010 + re_a100001)*ry2)*re2
                                        reI2 += (re_a010000 + (re_a011000 + re_a010100)*rx2 + (re_a010010 + re_a010001)*ry2)*re2
                                        reI2 += (re_a001000 + re_a001100*rx2 + (re_a001010 + re_a001001)*ry2)*rx2
                                        reI2 += (re_a000100 + (re_a000110 + re_a000101)*ry2)*rx2 + (re_a000010 + re_a000011*ry2)*ry2 + re_a000001*ry2 + re_a000000
                                        #imI2 = (im_a100000 + im_a110000*re2 + (im_a101000 + im_a100100)*rx2 + (im_a100010 + im_a100001)*ry2)*re2
                                        #imI2 += (im_a010000 + (im_a011000 + im_a010100)*rx2 + (im_a010010 + im_a010001)*ry2)*re2
                                        #imI2 += (im_a001000 + im_a001100*rx2 + (im_a001010 + im_a001001)*ry2)*rx2
                                        #imI2 += (im_a000100 + (im_a000110 + im_a000101)*ry2)*rx2 + (im_a000010 + im_a000011*ry2)*ry2 + im_a000001*ry2 + im_a000000
                                        #OC05052018: Note that this does not include all terms of multi-dim "bi-linear" interpolation!?

                                        resDegCoh[resOfst] = absMI/(sqrt(reI1*reI2) + absZerTolI)
                                    else:
                                        resDegCoh[resOfst] = 0
                                else:
                                    ie1_0 = ie
                                    ie2_0 = iep
                                    ie1_0_origPerE = ie1_0*origPerE
                                    ie2_0_origPerEp = ie2_0*origPerEp
                                    ie2_0_origPerEp_p_ie1_0_origPerE = ie2_0_origPerEp + ie1_0_origPerE
                                    ie1_0_origPerEp_p_ie1_0_origPerE = ie1_0*origPerEp + ie1_0_origPerE
                                    ie2_0_origPerEp_p_ie2_0_origPerE = ie2_0_origPerEp + ie2_0*origPerE

                                    ofstMI = ie2_0_origPerEp_p_ie1_0_origPerE + ix2_0*origPerXp + ix1_0*origPerX + iy2_0*origPerYp + iy1_0*origPerY
                                    ofstI1 = ie1_0_origPerEp_p_ie1_0_origPerE + ix1_0*origPerXp + ix1_0*origPerX + iy1_0*origPerYp + iy1_0*origPerY
                                    ofstI2 = ie2_0_origPerEp_p_ie2_0_origPerE + ix2_0*origPerXp + ix2_0*origPerX + iy2_0*origPerYp + iy2_0*origPerY

                                    reMI = self.arS[ofstMI]; imMI = self.arS[ofstMI + 1]
                                    absMI = sqrt(reMI*reMI + imMI*imMI)
                                    reI1 = self.arS[ofstI1]#; imI1 = self.arS[ofstI1 + 1]
                                    reI2 = self.arS[ofstI2]#; imI2 = self.arS[ofstI2 + 1]
                                    resDegCoh[resOfst] = absMI/(sqrt(reI1*reI2) + absZerTolI)

                                resOfst += 1
                                resEp += resEstep
                            resE += resEstep
                        resXp += resXstep
                    resX += resXstep
                resYp += resYstep
            resY += resYstep

        return resDegCoh

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

    def delE(self, _type=0, _treatEX=1, _treatEY=1): #OC151115
        """Delete Electric Field data
        :param _type: type of data to be deleted: 0- arEx, arEy, arExAux, arEyAux; 1- arEx, arEy only; 2- arExAux, arEyAux only
        :param _treatEX: switch specifying whether Ex data should be deleted or not (1 or 0)
        :param _treatEY: switch specifying whether Ey data should be deleted or not (1 or 0)
        """
        if _treatEX:
            if((_type == 0) or (_type == 1)):
                if(self.arEx is not None):
                    del self.arEx; self.arEx = None
            if((_type == 0) or (_type == 2)):
                if(hasattr(self, 'arExAux')):
                    if(self.arExAux is not None):
                        del self.arExAux; self.arExAux = None
        if _treatEY:
            if((_type == 0) or (_type == 1)):
                if(self.arEy is not None):
                    del self.arEy; self.arEy = None
            if((_type == 0) or (_type == 2)):
                if(hasattr(self, 'arEyAux')):
                    if(self.arEyAux is not None):
                        del self.arEyAux; self.arEyAux = None

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

    def calc_stokes(self, _stokes, _n_stokes_comp=4): #OC04052018
    #def calc_stokes(self, _stokes):
        """Calculate Stokes parameters from Electric Field"""
        if(_stokes.mutual <= 0):
            nTot = self.mesh.ne*self.mesh.nx*self.mesh.ny
            #if(type(_stokes).__name__ != 'SRWLStokes')):
            if(isinstance(_stokes, SRWLStokes) == False):
                raise Exception("Incorrect Stokes parameters object submitted") 
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
                        auxArEx[ofstR + 1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty

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
                        auxArEy[ofstR + 1] = a000 + (a100 + (a110 + a111*ty)*tx + a101*ty)*te + (a010 + a011*ty)*tx + a001*ty
                        
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
        
#****************************************************************************
class SRWLOpt(object):
    """Optical Element (base class)"""

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
    
    def __init__(self, _nZones=100, _rn=0.1e-03, _thick=10e-06, _delta1=1e-06, _atLen1=0.1, _delta2=0, _atLen2=1e-06, _x=0, _y=0):
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
    
    def __init__(self, _nx=1, _ny=1, _rx=1e-03, _ry=1e-03, _arTr=None, _extTr=0, _Fx=1e+23, _Fy=1e+23, _x=0, _y=0, _ne=1, _eStart=0, _eFin=0):
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
            self.allocate(_ne, _nx, _ny)

        #self.ne = _ne #number of transmission data points vs photon energy
        #self.nx = _nx #numbers of transmission data points in the horizontal and vertical directions
        #self.ny = _ny
        #self.eStart = _eStart #initial and final values of photon energy
        #self.eFin = _eFin
        #self.rx = _rx #ranges of horizontal and vertical coordinates [m] for which the transmission is defined
        #self.ry = _ry

        halfRangeX = 0.5*_rx;
        halfRangeY = 0.5*_ry;
        self.mesh = SRWLRadMesh(_eStart, _eFin, _ne, _x - halfRangeX, _x + halfRangeX, _nx, _y - halfRangeY, _y + halfRangeY, _ny)

        self.extTr = _extTr #0- transmission outside the grid/mesh is zero; 1- it is same as on boundary
        self.Fx = _Fx #estimated focal lengths [m]
        self.Fy = _Fy
        
        #self.x = _x #transverse coordinates of center [m]
        #self.y = _y
        #if _ne > 1: _Fx, _Fy should be arrays vs photon energy?

    def allocate(self, _ne, _nx, _ny):
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
        self.arTr = array('d', [0]*nTot)

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
    
        if(not(isinstance(_refl, list) or isinstance(_refl, array))):
            self.arRefl = array('d', [_refl]*nTot)
            for i in range(int(round(nTot/2))):
                i2 = i*2
                self.arRefl[i2] = _refl
                self.arRefl[i2 + 1] = 0
        else:
            self.arRefl = _refl

        #DEBUG
        #print(self.arRefl)
        
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
        self.Fx = 0 #i.e. focal lengths are not set
        self.Fy = 0

    def set_all(self,
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

        self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

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
    
    def __init__(self, _mirSub, _m=1, _grDen=100, _grDen1=0, _grDen2=0, _grDen3=0, _grDen4=0, _grAng=0):
        """
        :param _mirSub: SRWLOptMir (or derived) type object defining substrate of the grating
        :param _m: output (diffraction) order
        :param _grDen: groove density [lines/mm] (coefficient a0 in the polynomial groove density: a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4)
        :param _grDen1: groove density polynomial coefficient a1 [lines/mm^2]
        :param _grDen2: groove density polynomial coefficient a2 [lines/mm^3]
        :param _grDen3: groove density polynomial coefficient a3 [lines/mm^4]
        :param _grDen4: groove density polynomial coefficient a4 [lines/mm^5]
        :param _grAng: angle between the grove direction and the saggital direction of the substrate [rad] (by default, groves are made along saggital direction (_grAng=0))
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

class SRWLOptCryst(SRWLOpt):
    """Optical Element: Ideal Crystal"""

    #def _init_(self, _d_space, _psiOr, _psiOi, _psiHr, _psiHi, _psiHBr, _psiHBi, _H1, _H2, _H3, _Tc, _Tasym, _nx, _ny, _nz, _sx, _sy, _sz, _aChi, _aPsi, _aThe):
    #def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _h1, _h2, _h3, _tc, _ang_as, _nvx, _nvy, _nvz, _tvx, _tvy):
    #def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0):
    def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _uc=1):
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

        self.aux_energy = None  #MR01082016: renamed self.energy to self.aux_energy.
        self.aux_ang_dif_pl = None  #MR01082016: renamed self.ang_dif_pl to self.aux_ang_dif_pl. 
        
        self.uc = _uc #OC04092016

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

    def find_orient(self, _en, _ang_dif_pl=0):
        """Finds optimal crystal orientation in the input beam frame (i.e. surface normal and tangential vectors) and the orientation of the output beam frame (i.e. coordinates of the longitudinal and horizontal vectors in the input beam frame)
        :param _en: photon energy [eV]
        :param _ang_dif_pl: diffraction plane angle (0 corresponds to the vertical deflection; pi/2 to the horizontal deflection; any value in between is allowed)
        :return: list of two triplets of vectors:
                out[0] is the list of 3 base vectors [tangential, saggital, normal] defining the crystal orientation
                out[1] is the list of 3 base vectors [ex, ey, ez] defining orientation of the output beam frame
                the cartesian coordinates of all these vectors are given in the frame of the input beam
        """

        self.aux_energy = _en #MR01082016: renamed self.energy to self.aux_energy.
        self.aux_ang_dif_pl = _ang_dif_pl #MR01082016: renamed self.ang_dif_pl to self.aux_ang_dif_pl. 
   
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

        def prodV(_a, _b):
            return [_a[1]*_b[2] - _a[2]*_b[1], _a[2]*_b[0] - _a[0]*_b[2], _a[0]*_b[1] - _a[1]*_b[0]]
        
        def prodMV(_m, _v):
            return [_m[0][0]*_v[0] + _m[0][1]*_v[1] + _m[0][2]*_v[2],
                _m[1][0]*_v[0] + _m[1][1]*_v[1] + _m[1][2]*_v[2],
                _m[2][0]*_v[0] + _m[2][1]*_v[1] + _m[2][2]*_v[2]]
        
        def normV(_a):
            return sqrt(sum(n**2 for n in _a))

        #Crystal orientation vectors
        nv = [0, cos(tIn), -sin(tIn)]
        tv = [0, sin(tIn), cos(tIn)]
        sv = prodV(nv, tv)

        mc = [[sv[0], nv[0], tv[0]],
              [sv[1], nv[1], tv[1]],
              [sv[2], nv[2], tv[2]]]

        #Diffracted beam frame vectors
        z1c = [sv[2], sqrt(1. - sv[2]**2 - (tv[2] + wA*hv[2])**2), tv[2] + wA*hv[2]]
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
        if(self.uc == 2): #Bragg Transmission
            ex = [1, 0, 0]
            ey = [0, 1, 0]
            ez = [0, 0, 1]  
        #To check this and implement other self.uc cases!

        return [[tvNew, svNew, nvNew], [ex, ey, ez]]

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
    print('Optical Element Setup: CRL Focal Length:', foc_len, 'm')

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

    for iy in range(_mask_Ny):
        # Calculate the relative position in y.
        # NOTE: Use round to solve the precision issue!
        pitch_num_y = floor(round(y / _pitch_y, 9))
        y_rel = y - (pitch_num_y * _pitch_y) - _mask_y0
        if y_rel >= _pitch_y / 2:
            y_rel -= _pitch_y

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

                if (k2 * x_rel + (yCross2 - k2 * xCross2)) > y_rel > (k3 * x_rel + (yCross1 - k3 * xCross1)) \
                        and (k1 * x_rel + (yCross1 - k1 * xCross1)) > y_rel > (k4 * x_rel + (yCross2 - k4 * xCross2)) \
                        and not (abs(x - _mask_x0) < _pitch_x / 2 and abs(y - _mask_y0) < _pitch_y / 2) \
                        and abs(x) < grid_Rx / 2 and abs(y) < grid_Ry / 2:
                    inside_hole = True

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

    #OCTEST
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

    #OCTEST
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
        sizeLong = _height_prof_data[0][npData - 1] - _height_prof_data[0][1]
    else:
        npData = len(_ar_arg_long)
        sizeLong = _ar_arg_long[npData - 1] - _ar_arg_long[0]
        
    sizeLongProj = sizeLong*sinAngR

    if _ar_arg_tr is None:
        npDataTr = len(_height_prof_data) - 1
        sizeTr = _height_prof_data[npDataTr - 1][0] - _height_prof_data[1][0]
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

    optSlopeErr = SRWLOptT(nx, ny, sizeX, sizeY)
    
    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)

    #print(auxMesh.xStart, auxMesh.xFin, auxMesh.nx, xStep)
    #print(xStep, yStep)

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
def srwl_opt_setup_gen_transm(_func_path_dif, _delta, _atten_len, _rx, _ry, _xc=0, _yc=0, _ext_tr=0, _fx=0, _fy=0, _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
    """
    Setup Transmission type Optical Element similar to one simulating CRL, but with arbitrary optical path difference over hor. and vert. positions, defined by external function _func_path_dif(x, y)
    :param _func_path_dif: user-defined function of 2 variables specifying optical path difference over hor. and vert. positions (x, y)
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
    y = -0.5*_ry #CRL is always centered on the grid, however grid can be shifted
    for iy in range(_ny):
        x = -0.5*_rx
        for ix in range(_nx):
            pathDif = _func_path_dif(x, y)
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
            fm1 = _func_path_dif(_xc - dx, _yc)
            f0 = _func_path_dif(_xc, _yc)
            f1 = _func_path_dif(_xc + dx, _yc)
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
            fm1 = _func_path_dif(_xc, _yc - dy)
            f0 = _func_path_dif(_xc, _yc)
            f1 = _func_path_dif(_xc, _yc + dy)
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

        resInt = SRWLStokes(1, 'f', self.eStartEff, self.eFinEff, 1, self.xStart, self.xFin, self.nx, self.yStart, self.yFin, self.ny)

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
                
            sktIn = SRWLStokes(arI, 'f', meshIn.eStart, meshIn.eFin, meshIn.ne, meshIn.xStart, meshIn.xFin, meshIn.nx, meshIn.yStart, meshIn.yFin, meshIn.ny)
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
                    effMult = interp_1d(ePh, self.eStartEff, self.eStepEff, self.neEff, self.specEff, _ord=2)
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

    if(_mutual != 0):
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
        #strOut += '#' + repr(_n_stokes) + ' #Number of components\n'
        nComp = _n_stokes
    nRadPt = _mesh.ne*_mesh.nx*_mesh.ny
    if(_mutual > 0): nRadPt *= nRadPt
    
    nVal = nRadPt*nComp #_mesh.ne*_mesh.nx*_mesh.ny*nComp
    if(_cmplx != 0): nVal *= 2 #OC06052018

    for i in range(nVal): #write all data into one column using "C-alignment" as a "flat" 1D array
        f.write(' ' + repr(_ar_intens[i]) + '\n')
        #strOut += ' ' + repr(_ar_intens[i]) + '\n'

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
        text_to_save = '[{}]: Calculation on {} cores with averaging of {} particles/iteration.'.format(
            timestamp,
            cores,
            particles_per_iteration,
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
def srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, num_of_proc, num_part_per_proc, num_sent_per_proc, num_part_avg_proc): #MR20012017  
    """Initialize parameters for the SRW status files and generate the files.  
    :param rank: rank of the process.  
    :param num_of_proc: total number of processes.  
    :param num_part_per_proc: number of electrons treated by each worker process.  
    :param num_sent_per_proc: number of sending acts made by each worker process.  
    :param num_part_avg_proc: number of macro-electrons to be used in calculation by each worker.  
    """  
    log_dir = os.path.abspath('__srwl_logs__')  
    if rank == 0:  
        try:  
            os.mkdir(log_dir)  
        except OSError as exc:  
            if exc.errno == errno.EEXIST and os.path.isdir(log_dir):
                pass  
            else:
                raise  
    timestamp = '{:%Y-%m-%d_%H-%M-%S}'.format(datetime.datetime.now())  
    log_file = 'srwl_stat_wfr_emit_prop_multi_e_{}'.format(timestamp)  
    log_path = os.path.join(log_dir, log_file)  
    if num_of_proc <= 1:  
        total_num_of_particles = num_part_per_proc  
    else:  
        total_num_of_particles = num_sent_per_proc * (num_of_proc - 1) * num_part_avg_proc  
    if rank == 0:  
        srwl_uti_save_stat_wfr_emit_prop_multi_e(0, total_num_of_particles, filename=log_path, cores=num_of_proc, particles_per_iteration=num_part_avg_proc)  
    return log_path, total_num_of_particles

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
def srwl_uti_array_alloc(_type, _n):
    nPartMax = 10000000 #to tune
    if(_n <= nPartMax): return array(_type, [0]*_n)
        #resAr = array(_type, [0]*_n)
        #print('Array requested:', _n, 'Allocated:', len(resAr))
        #return resAr

    nEqualParts = int(_n/nPartMax)
    nResid = int(_n - nEqualParts*nPartMax)
    resAr = array(_type, [0]*nPartMax)
    if(nEqualParts > 1):
        auxAr = deepcopy(resAr)
        for i in range(nEqualParts - 1): resAr.extend(auxAr)
    if(nResid > 0):
        auxAr = array(_type, [0]*nResid)
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
    
    for iz in range(0, nzp1):

        if(iz > 0):
            print('Propagation (step # ' + repr(iz) + ') ... ', end='')
            t0 = time.time();
            srwl.PropagElecField(_wfr, cntDrift)
            print('completed (lasted', round(time.time() - t0, 6), 's)')

        curMesh = _wfr.mesh
        xStepCurMesh = (curMesh.xFin - curMesh.xStart)/(curMesh.nx - 1)
        yStepCurMesh = (curMesh.yFin - curMesh.yStart)/(curMesh.ny - 1)
        curIntVsX = array('f', [0]*_wfr.mesh.nx)
        curIntVsY = array('f', [0]*_wfr.mesh.ny)
        
        srwl.CalcIntFromElecField(curIntVsX, _wfr, _pol, _type, 1, ec, _xc, _yc) #int. vs X

        x = xStart
        for ix in range(_nx): #interpolation
            resIntVsZX[iz + ix*nzp1] = uti_math.interp_1d(x, curMesh.xStart, xStepCurMesh, curMesh.nx, curIntVsX, _ord_interp)
            x += xStep
        
        srwl.CalcIntFromElecField(curIntVsY, _wfr, _pol, _type, 2, ec, _xc, _yc) #int. vs Y

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

    return resIntVsZX, resIntVsZY, resIntVsZXY, resMesh

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
    wfr.unitElFld = 1 #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time) ?
    #wfr.arElecPropMatr = array('d', [0] * 20) #effective 1st order "propagation matrix" for electron beam parameters
    #wfr.arMomX = array('d', [0] * 11 * _ne) #statistical moments (of Wigner distribution); to check the exact number of moments required
    #wfr.arMomY = array('d', [0] * 11 * _ne)
    #wfr.arWfrAuxData = array('d', [0] * 30) #array of auxiliary wavefront data

    return wfr

#**********************Main Partially-Coherent Emission and Propagaiton simulation function
def srwl_wfr_emit_prop_multi_e(_e_beam, _mag, _mesh, _sr_meth, _sr_rel_prec, _n_part_tot, _n_part_avg_proc=1, _n_save_per=100,
                               _file_path=None, _sr_samp_fact=-1, _opt_bl=None, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0,
                               #_rand_meth=1, _tryToUseMPI=True):
                               #_rand_meth=1, _tryToUseMPI=True, _w_wr=0.): #OC26032016 (added _w_wr)
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0.): #OC07092016 (renamed _wr)
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _det=None): #OC06122016
                               #_rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None): #OC05012017
                               _rand_meth=1, _tryToUseMPI=True, _wr=0., _wre=0., _det=None, _me_approx=0): #OC05042017
    """
    Calculate Stokes Parameters of Emitted (and Propagated, if beamline is defined) Partially-Coherent SR
    :param _e_beam: Finite-Emittance e-beam (SRWLPartBeam type)
    :param _mag: Magnetic Field container (magFldCnt type) or Coherent Gaussian beam (SRWLGsnBm type) or Point Source (SRWLPtSrc type)
    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
    :param _sr_meth: SR Electric Field calculation method to be used: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
    :param _sr_rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, the smaller the more accurate)
    :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
    :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
    :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
    :param _file_path: path to file for saving intermediate average Stokes data by master process
    :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
    :param _opt_bl: optical beamline (container) to propagate the radiation through (SRWLOptC type)
    :param _pres_ang: switch specifying presentation of the resulting Stokes parameters: coordinate (0) or angular (1)
    :param _char: radiation characteristic to calculate:
        0- Total Intensity, i.e. Flux per Unit Surface Area (s0);
        1- Four Stokes components of Flux per Unit Surface Area;
        2- Mutual Intensity Cut vs X;
        3- Mutual Intensity Cut vs Y;
        4- Mutual Intensity Cuts and Degree of Coherence vs X & Y;
        10- Flux
        20- Electric Field (sum of fields from all macro-electrons, assuming CSR)
        40- Total Intensity, i.e. Flux per Unit Surface Area (s0), Mutual Intensity Cuts and Degree of Coherence vs X & Y;
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
    """

    doMutual = 0 #OC30052017
    if((_char >= 2) and (_char <= 4)): doMutual = 1

    if((_det is not None) and (doMutual > 0)): raise Exception("Detector processing is not supported for mutual intensity") #OC30052017

    import time #DEBUG
    #print('_mesh.xStart=', _mesh.xStart, '_mesh.xFin=', _mesh.xFin, '_mesh.yStart=', _mesh.yStart, '_mesh.yFin=', _mesh.yFin) #DEBUG
    #print('_sr_samp_fact:', _sr_samp_fact) #DEBUG
    #for i in range(len(_opt_bl.arProp)): #DEBUG
    #    print(i, _opt_bl.arProp[i]) #DEBUG

    #DEBUG
    #self.arOpt = []
    #self.arProp = []

    nProc = 1
    rank = 1
    MPI = None
    comMPI = None

    if(_tryToUseMPI):
        try:
            ##DEBUG
            ##resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], -1) #MPI module load
            #resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], 0) #MPI module load
            ##print('__import__ passed')
            #MPI = resImpMPI4Py.MPI

            from mpi4py import MPI #OC091014
        
            comMPI = MPI.COMM_WORLD
            rank = comMPI.Get_rank()
            nProc = comMPI.Get_size()

        except:
            print('Calculation will be sequential (non-parallel), because "mpi4py" module can not be loaded')

    #print('DEBUG:', MPI)
    #print('DEBUG: rank, nProc:', rank, nProc)

    #if(nProc <= 1): #OC050214
    #    _n_part_avg_proc = _n_part_tot

    #OC30052017 (commented-out)
    #wfr = SRWLWfr() #Wavefronts to be used in each process
    #wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    #wfr.mesh.set_from_other(_mesh)    
    #wfr.partBeam = deepcopy(_e_beam)

    #arPrecParSR = [_sr_meth, _sr_rel_prec, 0, 0, 50000, 0, _sr_samp_fact] #to add npTraj, useTermin ([4], [5]) terms as input parameters
    arPrecParSR = [_sr_meth, _sr_rel_prec, 0, 0, 50000, 1, _sr_samp_fact] #to add npTraj, useTermin ([4], [5]) terms as input parameters

    meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) #OC30052017 (uncommented) #to ensure correct final mesh if _opt_bl is None
    #meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) if(_det is None) else _det.get_mesh() #OC06122016
    #Note: the case ((_det is not None) and (doMutual > 0)) is not supported

    file_path1 = copy(_file_path) #OC30052017
    file_path2 = None
    file_path_deg_coh1 = None #07052018
    file_path_deg_coh2 = None
    meshRes2 = None #OC30052017 
    #if(_char == 4): #Mutual Intensity Cuts vs X & Y are required
    if((_char == 4) or (_char == 40)): #OC02052018 #Mutual Intensity Cuts vs X & Y are required
        #meshRes2 = copy(meshRes)
        if(_char == 4): meshRes2 = copy(meshRes) #OC02052018
        file_path1 += '.1'
        file_path2 = copy(_file_path) + '.2'
        file_path_deg_coh1 = copy(_file_path) + '.dc.1'
        file_path_deg_coh2 = copy(_file_path) + '.dc.2'

    wfr = SRWLWfr() #Wavefronts to be used in each process
    wfr2 = None #OC30052017
    if((_opt_bl is None) and (doMutual > 0)): #OC30052017
        if(_char == 2): #Cut vs X
            meshRes.ny = 1
            meshRes.yStart = _y0
            meshRes.yFin = _y0
        elif(_char == 3): #Cut vs Y
            meshRes.nx = 1
            meshRes.xStart = _x0
            meshRes.xFin = _x0
        elif(_char == 4): #Cuts of Mutual Intensity vs X & Y
            meshRes.ny = 1
            meshRes.yStart = _y0
            meshRes.yFin = _y0
            meshRes2.nx = 1
            meshRes2.xStart = _x0
            meshRes2.xFin = _x0
            wfr2 = SRWLWfr()

    #OC30052017
    wfr.allocate(meshRes.ne, meshRes.nx, meshRes.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    wfr.mesh.set_from_other(meshRes)    
    wfr.partBeam = deepcopy(_e_beam)
    if(wfr2 is not None):
        wfr2.allocate(meshRes2.ne, meshRes2.nx, meshRes2.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr2.mesh.set_from_other(meshRes2)    
        wfr2.partBeam = deepcopy(_e_beam)

    if(_det is not None): meshRes = _det.get_mesh() #OC06122016

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
    if((_char == 10) and (_mesh.nx == 1) and (_mesh.ny == 1)):
        calcSpecFluxSrc = True
        ePhIntegMult *= 1.e+06*(_mesh.xFin - _mesh.xStart)*(_mesh.yFin - _mesh.yStart) #to obtain Flux from Intensity (Flux/mm^2)

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

    #_sr_rel_prec = int(_sr_rel_prec)
    
    _n_part_tot = int(_n_part_tot)
    _n_part_avg_proc = int(_n_part_avg_proc)
    if(_n_part_avg_proc <= 0): _n_part_avg_proc = 1
    _n_save_per = int(_n_save_per)

    nPartPerProc = _n_part_tot
    nSentPerProc = 0
    
    if(nProc <= 1):
        _n_part_avg_proc = _n_part_tot
    else: #OC050214: adjustment of all numbers of points, to make sure that sending and receiving are consistent
 
        nPartPerProc = int(round(_n_part_tot/(nProc - 1)))
        nSentPerProc = int(round(nPartPerProc/_n_part_avg_proc)) #Number of sending acts made by each worker process

        if(nSentPerProc <= 0): #OC160116
            nSentPerProc = 1
            _n_part_avg_proc = nPartPerProc
        
        nPartPerProc = _n_part_avg_proc*nSentPerProc #Number of electrons treated by each worker process

    #print('DEBUG MESSAGE: rank:', rank,': nPartPerProc=', nPartPerProc, 'nSentPerProc=', nSentPerProc, '_n_part_avg_proc=', _n_part_avg_proc)

    log_path, total_num_of_particles = srwl_uti_save_stat_wfr_emit_prop_multi_e_init(rank, nProc, nPartPerProc, nSentPerProc, _n_part_avg_proc) #MR20012017

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
    resStokes2 = None #OC30052017
    workStokes2 = None
    arAuxResSt12 = None #OC31052017
    arAuxWorkSt12 = None #OC31052017

    resStokes3 = None #OC03052018
    workStokes3 = None
    arAuxResSt123 = None #OC03052018
    arAuxWorkSt123 = None #OC03052018

    iAvgProc = 0
    iSave = 0

    #doMutual = 0
    #if((_char >= 2) and (_char <= 4)): doMutual = 1

    #depTypeME_Approx = 3 #vs x&y (horizontal and vertical positions or angles); to be used only at _me_approx > 0
    #OC15092017
    depTypeInt = 3 #vs x&y (horizontal and vertical positions or angles); to be used only at _me_approx > 0
    #phEnME_Approx = _mesh.eStart
    phEnInt = _mesh.eStart #OC15092017
    #if(_me_approx == 1): #OC15092017
    #if((_mesh.nx <= 1) or (_mesh.ny <= 1)):
    if((_me_approx == 1) and ((_mesh.nx <= 1) or (_mesh.ny <= 1))): #OC01102017
        raise Exception("This calculation method requires more than one observation point in the horizontal and vertical directions")

    if(_mesh.ne > 1):
        #OC15092017
        depTypeInt = 6 #vs e&x&y (photon energy or time, horizontal and vertical positions or angles)
        phEnInt = 0.5*(_mesh.eStart + _mesh.eFin)
        #depTypeME_Approx = 6 #vs e&x&y (photon energy or time, horizontal and vertical positions or angles)
        #phEnME_Approx = 0.5*(_mesh.eStart + _mesh.eFin)

    resEntityName = 'Intensity' #OC26042016
    resEntityUnits = 'ph/s/.1%bw/mm^2'
    resLabelsToSave = ['Photon Energy', 'Horizontal Position', 'Vertical Position', resEntityName] #OC26042016
    resUnitsToSave = ['eV', 'm', 'm', resEntityUnits] #OC26042016
    
    if(_pres_ang > 0): #OC20112017
        resEntityName = 'Ang. Intensity Distr.'
        resEntityUnits = 'ph/s/.1%bw/mrad^2'
        resLabelsToSave = ['Photon Energy', 'Horizontal Angle', 'Vertical Angle', resEntityName]
        resUnitsToSave = ['eV', 'rad', 'rad', resEntityUnits]

    if(calcSpecFluxSrc == True):
        resEntityName = 'Flux'
        resEntityUnits = 'ph/s/.1%bw'
        resLabelsToSave[3] = resEntityName #OC27032018
        resUnitsToSave[3] = resEntityUnits #OC27032018
    
    #resLabelsToSaveMutualHorCut = [resLabelsToSave[0], resLabelsToSave[1], 'Conj. ' + resLabelsToSave[1], 'Mutual ' + resLabelsToSave[3]] #OC03052018
    #resLabelsToSaveMutualVerCut = [resLabelsToSave[0], resLabelsToSave[2], 'Conj. ' + resLabelsToSave[2], 'Mutual ' + resLabelsToSave[3]]

    #if(((rank == 0) or (nProc == 1)) and (_opt_bl != None)): #calculate once the central wavefront in the master process (this has to be done only if propagation is required)
    if(((rank == 0) or (nProc == 1)) and (_det is None)): #12/01/2017 

        if(useGsnBmSrc):
            srwl.CalcElecFieldGaussian(wfr, _mag, arPrecParSR)
            if(wfr2 is not None): srwl.CalcElecFieldGaussian(wfr2, _mag, arPrecParSR) #OC30052017
            #print('DEBUG: Commented-out: CalcElecFieldGaussian')
        elif(usePtSrc): #OC16102017
            srwl.CalcElecFieldPointSrc(wfr, _mag, arPrecParSR)
            if(wfr2 is not None): srwl.CalcElecFieldPointSrc(wfr2, _mag, arPrecParSR)
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
            srwl.PropagElecField(wfr, _opt_bl)

        #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
        #meshRes.set_from_other(wfr.mesh) #DEBUG
        #resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual) #DEBUG
        #wfr.calc_stokes(resStokes) #DEBUG
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual) #DEBUG

        #print('DEBUG: Commented-out: PropagElecField')
        #print('DEBUG MESSAGE: Central Wavefront propagated')
        if(_pres_ang > 0):
            wfr.unitElFldAng = 1 #OC20112017 (to ensure result in [ph/s/.1%bw/mrad^2])
            srwl.SetRepresElecField(wfr, 'a')
            if(wfr2 is not None):
                wfr2.unitElFldAng = 1 #OC20112017
                srwl.SetRepresElecField(wfr2, 'a') #OC30052017
            #print('DEBUG: Commented-out: SetRepresElecField')

        #meshRes.set_from_other(wfr.mesh)
        #if(_det is None): meshRes.set_from_other(wfr.mesh) #OC06122016
        #else: meshRes = _det.get_mesh() #??
        
        meshRes.set_from_other(wfr.mesh) #OC12012016

        #if(wfr2 is not None): meshRes2.set_from_other(wfr2.mesh) #OC30052017

        if(doMutual > 0):
            if(_char == 2): #Cut vs X
                meshRes.ny = 1
                meshRes.yStart = _y0
                meshRes.yFin = _y0
            elif(_char == 3): #Cut vs Y
                meshRes.nx = 1
                meshRes.xStart = _x0
                meshRes.xFin = _x0
            elif(_char == 4): #Cuts of Mutual Intensity vs X & Y
                meshRes.ny = 1
                meshRes.yStart = _y0
                meshRes.yFin = _y0
                if(wfr2 is None): meshRes2.set_from_other(wfr.mesh) #OC31052017
                else: meshRes2.set_from_other(wfr2.mesh) #OC31052017
                meshRes2.nx = 1
                meshRes2.xStart = _x0
                meshRes2.xFin = _x0

        if(nProc > 1): #send resulting mesh to all workers
            #comMPI.send(wfr.mesh, dest=)
            arMesh = array('f', [meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny])
            
            #arMesh2 = None #OC30052017
            if(_char == 4): #Cuts of Mutual Intensity vs X & Y
                arMesh2 = array('f', [meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny])
                for ii in range(len(arMesh2)): #OC31052017
                    arMesh.append(arMesh2[ii])

            #comMPI.Bcast([arMesh, MPI.FLOAT], root=MPI.ROOT)
            #comMPI.Bcast([arMesh, MPI.FLOAT])

            #print('DEBUG MESSAGE: Rank0 is about to broadcast mesh of Propagated central wavefront')
            for iRank in range(nProc - 1):
                dst = iRank + 1
                #print("msg %d: sending data from %d to %d" % (iRank, rank, dst)) #an he
                comMPI.Send([arMesh, MPI.FLOAT], dest=dst)

                #OC30052017
                #if(_char == 4): #Cuts of Mutual Intensity vs X & Y
                #    comMPI.Send([arMesh2, MPI.FLOAT], dest=dst)

            #print('DEBUG MESSAGE: Mesh of Propagated central wavefront broadcasted')

        #DEBUG
        #print('meshRes: ne=', meshRes.ne, 'eStart=', meshRes.eStart, 'eFin=', meshRes.eFin)
        #print('meshRes: nx=', meshRes.nx, 'xStart=', meshRes.xStart, 'xFin=', meshRes.xFin)
        #print('meshRes: ny=', meshRes.ny, 'yStart=', meshRes.yStart, 'yFin=', meshRes.yFin)
        #END DEBUG

        resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
        #wfr.calc_stokes(resStokes) #OC190414: don't take into account the first "central" beam!
        #workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
        #OC06042017 (commented-out workStokes = ...)

        if(_char == 4): #OC31052017 #Cuts of Mutual Intensity vs X & Y
            resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
            lenArSt12 = len(resStokes.arS) + len(resStokes2.arS)
            arAuxResSt12 = array('f', [0]*lenArSt12)

        if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
            resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
            resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
            lenArSt123 = len(resStokes.arS) + len(resStokes2.arS) + len(resStokes3.arS)
            arAuxResSt123 = array('f', [0]*lenArSt123)

        #iAvgProc += 1 #OC190414 (commented-out)
        #iSave += 1

    #slaves = [] #an he
    #print('DEBUG MESSAGE: rank=', rank)
    numComp = 1
    if(_char == 20): numComp = 4 #OC16012017

    if(_opt_bl is None): arPrecParSR[6] = 0 #Ensure non-automatic choice of numbers of points if there is no beamline

    if((rank > 0) or (nProc == 1)):

        #if((nProc > 1) and (_opt_bl != None)): #receive mesh for the resulting wavefront from the master
        #if((nProc > 1) and (_opt_bl is not None) and (_det is None)): #OC12012017 #receive mesh for the resulting wavefront from the master
        if((nProc > 1) and (_det is None)): #OC16012017 #receive mesh for the resulting wavefront from the master
            #arMesh = array('f', [0]*9)
            nNumToRecv = 9
            if(_char == 4): nNumToRecv = 18 #OC31052017 #Cuts of Mutual Intensity vs X & Y
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

            #arMesh2 = None #OC30052017
            if(_char == 4): #Cuts of Mutual Intensity vs X & Y
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

            #sys.exit(0)

        nRadPt = meshRes.ne*meshRes.nx*meshRes.ny
        if(doMutual > 0): nRadPt *= nRadPt
        nStPt = nRadPt*4

        nRadPt2 = None #OC30052017
        nStPt2 = None
        if(_char == 4): #Cuts of Mutual Intensity vs X & Y
            nRadPt2 = meshRes2.ne*meshRes2.nx*meshRes2.ny
            if(doMutual > 0): nRadPt2 *= nRadPt2 #OC03052018
            nStPt2 = nRadPt2*4

        nRadPt3 = None #OC03052018
        nStPt3 = None
        if(_char == 40): #Intensity and Cuts of Mutual Intensity vs X & Y
            nRadPt2 = meshRes.ne*meshRes.nx
            nRadPt2 *= nRadPt2
            nStPt2 = nRadPt2*4
            nRadPt3 = meshRes.ne*meshRes.ny
            nRadPt3 *= nRadPt3
            nStPt3 = nRadPt3*4
        
        randAr = array('d', [0]*6) #for random Gaussian numbers

        #random.seed(rank)
        random.seed(rank*123)
        newSeed = random.randint(0, 1000000)
        random.seed(newSeed)
        
        iAuxSendCount = 0 #for debug
        
        for i in range(nPartPerProc): #loop over macro-electrons

            if(_me_approx == 0): #OC05042017 #General method
                
                if(_rand_meth == 1):
                    for ir in range(5): #to expend to 6D eventually
                        randAr[ir] = random.gauss(0, 1)
                elif(_rand_meth == 2):
                    if(nProc > 1):
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
                #    randAr = array('d', [0,0,0,2,0])
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

            elif(_me_approx == 1): #OC05042017 #Numerical integration only over electron energy

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
                else: print('i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
                
                if(_e_ph_integ == 1): print('Eph=', wfr.mesh.eStart)

            if(calcSpecFluxSrc): #consider taking into account _rand_meth != 1 here
                xObs = random.uniform(_mesh.xStart, _mesh.xFin)
                wfr.mesh.xStart = xObs
                wfr.mesh.xFin = xObs
                yObs = random.uniform(_mesh.yStart, _mesh.yFin)
                wfr.mesh.yStart = yObs
                wfr.mesh.yFin = yObs
                #print('xObs=', xObs, 'yObs=', yObs)

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
                else:

                    #print('Single-electron SR calculatiton ... ', end='') #DEBUG
                    #t0 = time.time(); #DEBUG
                    #print('arPrecParSR[6]=', arPrecParSR[6])

                    #DEBUG
                    #print('Numbers of points (after re-setting):') #DEBUG
                    #print('ne=', wfr.mesh.ne, 'eStart=', wfr.mesh.eStart, 'eFin=', wfr.mesh.eFin)
                    #print('nx=', wfr.mesh.nx, 'xStart=', wfr.mesh.xStart, 'xFin=', wfr.mesh.xFin)
                    #print('ny=', wfr.mesh.ny, 'yStart=', wfr.mesh.yStart, 'yFin=', wfr.mesh.yFin)
                    #print('zStart=', wfr.mesh.zStart)
                    #END DEBUG

                    #DEBUG
                    #if(i == 48):
                    #    print('Eel=', wfr.partBeam.partStatMom1.get_E(), ' GeV')
                    #    print('xe=', wfr.partBeam.partStatMom1.x, ' xpe=', wfr.partBeam.partStatMom1.xp, ' ye=', wfr.partBeam.partStatMom1.y, ' ype=', wfr.partBeam.partStatMom1.yp)
                    #    print('nx=', wfr.mesh.nx, ' xStart=', wfr.mesh.xStart, ' xFin=', wfr.mesh.xFin)
                    #    print('ny=', wfr.mesh.ny, ' yStart=', wfr.mesh.yStart, ' yFin=', wfr.mesh.yFin)
                    #END DEBUG

                    srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR) #calculate Electric Field emitted by current electron
                    if(wfr2 is not None): srwl.CalcElecFieldSR(wfr2, 0, _mag, arPrecParSR) #OC30052017

                    #DEBUG
                    #if(i == 48):
                    #    arI1 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
                    #    srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
                    #    srwl_uti_save_intens_ascii(arI1, wfr.mesh, os.path.join(os.getcwd(), 'data_CDI', 'debug_int_se.dat'))
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

                    srwl.PropagElecField(wfr, _opt_bl) #propagate Electric Field emitted by the electron

                    #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                    #print('DEBUG: Commented-out: PropagElecField')

                    #DEBUG
                    #if(i == 48):
                    #    arI1 = array('f', [0]*wfr.mesh.nx*wfr.mesh.ny) #"flat" 2D array to take intensity data
                    #    srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, wfr.mesh.eStart, 0, 0)
                    #    srwl_uti_save_intens_ascii(arI1, wfr.mesh, os.path.join(os.getcwd(), 'data_CDI', 'debug_int_pr_se.dat'))
                    #    sys.exit()
                    #END DEBUG

                if(_pres_ang > 0):
                    wfr.unitElFldAng = 1 #OC20112017 (to have result in [ph/s/.1%bw/mrad^2] vs [rad])
                    srwl.SetRepresElecField(wfr, 'a')
                    if(wfr2 is not None):
                        wfr2.unitElFldAng = 1 #OC20112017
                        srwl.SetRepresElecField(wfr2, 'a')

                    #print('DEBUG: Commented-out: SetRepresElecField')

            except:
                traceback.print_exc()

            meshWork = deepcopy(wfr.mesh)
            meshWork2 = None #OC30052017
            if(doMutual > 0):
                if(_char == 2):
                    meshWork.ny = 1
                    meshWork.yStart = _y0
                    meshWork.yFin = _y0
                elif(_char == 3):
                    meshWork.nx = 1
                    meshWork.xStart = _x0
                    meshWork.xFin = _x0
                elif(_char == 4): #OC30052017 #Cuts of Mutual Intensity vs X & Y
                    meshWork.ny = 1
                    meshWork.yStart = _y0
                    meshWork.yFin = _y0
                    meshWork2 = deepcopy(wfr.mesh)
                    meshWork2.nx = 1
                    meshWork2.xStart = _x0
                    meshWork2.xFin = _x0

            #OC06042017 (commented-out the above, entered workStokes = ... below)
            workStokes = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, meshWork.yStart, meshWork.yFin, meshWork.ny, doMutual)
            
            if(_char == 4): #Cuts of Mutual Intensity vs X & Y
                workStokes2 = SRWLStokes(1, 'f', meshWork2.eStart, meshWork2.eFin, meshWork2.ne, meshWork2.xStart, meshWork2.xFin, meshWork2.nx, meshWork2.yStart, meshWork2.yFin, meshWork2.ny, doMutual)
                     
            if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
                workStokes2 = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, _y0, _y0, 1, _mutual=1)
                workStokes3 = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, _x0, _x0, 1, meshWork.yStart, meshWork.yFin, meshWork.ny, _mutual=1)

            if(_me_approx == 0): #OC05042017 #General case of numerical integration over 5D phase space of electron beam
                
                if(_char == 20): 
                    wfr.copy_comp(workStokes) #OC15012017: copy electric field components to Stokes structure
                elif(_char == 0): #OC15092017
                    srwl.CalcIntFromElecField(workStokes.arS, wfr, 6, 0, depTypeInt, phEnInt, 0., 0.)
                else: 
                    wfr.calc_stokes(workStokes, _n_stokes_comp=numComp) #calculate Stokes parameters from Electric Field
                    if(workStokes2 is not None): #OC30052017
                        if(wfr2 is None): wfr.calc_stokes(workStokes2, _n_stokes_comp=numComp)
                        else: wfr2.calc_stokes(workStokes2, _n_stokes_comp=numComp)

                    if(workStokes3 is not None): #OC03052018
                        wfr.calc_stokes(workStokes3, _n_stokes_comp=numComp)

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

                #elif(_char == 1): #Four Stokes components
                    #To implement extraction of Stokes components, with and without convolution, in CalcIntFromElecField

            #DEBUG
            #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, 1)
            #END DEBUG

            if(resStokes is None):
                resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
                #DEBUG
                #print('resStokes #2: ne=', resStokes.mesh.ne, 'eStart=', resStokes.mesh.eStart, 'eFin=', resStokes.mesh.eFin)
                #END DEBUG

            if(_char == 4): #OC31052017 #Cuts of Mutual Intensity vs X & Y
                if(resStokes2 is None):
                    resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
                if(arAuxResSt12 is None):
                    lenArSt12 = len(resStokes.arS) + len(resStokes2.arS)
                    arAuxResSt12 = array('f', [0]*lenArSt12)

            if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
                if(resStokes2 is None):
                    resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
                if(resStokes3 is None):
                    resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
                if(arAuxResSt123 is None):
                    lenArSt123 = len(resStokes.arS) + len(resStokes2.arS) + len(resStokes3.arS)
                    arAuxResSt123 = array('f', [0]*lenArSt123)

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
                
                #resStokes.avg_update_same_mesh(workStokes, iAvgProc, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
                resStokes.avg_update_same_mesh(workStokes, iAvgProc, numComp, ePhIntegMult) #OC16012017 #to treat all Stokes components / Polarization in the future

                if((resStokes2 is not None) and (workStokes2 is not None)): #OC30052017
                    resStokes2.avg_update_same_mesh(workStokes2, iAvgProc, numComp, ePhIntegMult)

                if((resStokes3 is not None) and (workStokes3 is not None)): #OC03052018
                    resStokes3.avg_update_same_mesh(workStokes3, iAvgProc, numComp, ePhIntegMult)

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

                if(_char == 40): #OC03052018
                    resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult)

                    if((resStokes2 is not None) and (workStokes2 is not None)): 
                        resStokes2.avg_update_interp_mutual(workStokes2, iAvgProc, 1, ePhIntegMult)

                    if((resStokes3 is not None) and (workStokes3 is not None)):
                        resStokes3.avg_update_interp_mutual(workStokes3, iAvgProc, 1, ePhIntegMult)

                else:
                    if(doMutual == 0): resStokes.avg_update_interp(workStokes, iAvgProc, 1, numComp, ePhIntegMult) #OC16012017 #to treat all Stokes components / Polarization in the future
                    else:
                        resStokes.avg_update_interp_mutual(workStokes, iAvgProc, 1, ePhIntegMult)

                        if((resStokes2 is not None) and (workStokes2 is not None)): #OC30052017
                            resStokes2.avg_update_interp_mutual(workStokes2, iAvgProc, 1, ePhIntegMult)

                #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                #print('DEBUG MESSAGE: Finished interpolation of current wavefront on resulting mesh')

            iAvgProc += 1
            if(iAvgProc >= _n_part_avg_proc):
                if(nProc > 1):
                    #sys.exit(0)
                    #print("sending data from %d to 0" % rank) #an he
                    #DEBUG
                    #srwl_uti_save_intens_ascii(resStokes.arS, resStokes.mesh, _file_path, 1)
                    #END DEBUG

                    #DEBUG
                    #srwl_uti_save_text("Preparing to sending # " + str(iAuxSendCount + 1), _file_path + "." + str(rank) + "bs.dbg")
                    #END DEBUG

                    #comMPI.Send([resStokes.arS, MPI.FLOAT], dest=0)
                    if(resStokes2 is None): comMPI.Send([resStokes.arS, MPI.FLOAT], dest=0)
                    #else: #OC31052017
                    elif(resStokes3 is None): #OC03052018
                        lenArSt1 = len(resStokes.arS)
                        lenArSt2 = len(resStokes2.arS)
                        for i1 in range(lenArSt1): arAuxResSt12[i1] = resStokes.arS[i1]
                        for i2 in range(lenArSt2): arAuxResSt12[i2 + lenArSt1] = resStokes2.arS[i2]
                        comMPI.Send([arAuxResSt12, MPI.FLOAT], dest=0)
                    else: #OC03052018
                        lenArSt1 = len(resStokes.arS)
                        lenArSt2 = len(resStokes2.arS)
                        lenArSt3 = len(resStokes3.arS)
                        for i1 in range(lenArSt1): arAuxResSt123[i1] = resStokes.arS[i1]
                        for i2 in range(lenArSt2): arAuxResSt123[i2 + lenArSt1] = resStokes2.arS[i2]
                        lenArSt12 = lenArSt1 + lenArSt2
                        for i3 in range(lenArSt3): arAuxResSt123[i3 + lenArSt12] = resStokes3.arS[i3]
                        comMPI.Send([arAuxResSt123, MPI.FLOAT], dest=0)

                    #if(resStokes2 is not None): comMPI.Send([resStokes2.arS, MPI.FLOAT], dest=0) #OC30052017

                    iAuxSendCount += 1 #for debug

                    #DEBUG
                    #srwl_uti_save_text("Sent # " + str(iAuxSendCount), _file_path + "." + str(rank) + "es.dbg")
                    #END DEBUG

                    for ir in range(nStPt):
                        resStokes.arS[ir] = 0

                    if(resStokes2 is not None): #OC30052017
                        for ir in range(nStPt2):
                            resStokes2.arS[ir] = 0

                    if(resStokes3 is not None): #OC03052018
                        for ir in range(nStPt3):
                            resStokes3.arS[ir] = 0

                    #DEBUG
                    #srwl_uti_save_intens_ascii(resStokes.arS, resStokes.mesh, _file_path, 1)
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

                    if(_char == 40): #OC03052018
                        srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0)
                        if(resStokes2 is not None):
                            #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                            srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                            srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC06052018
                        if(resStokes3 is not None):
                            #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                            srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                            srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC06052018
                    elif(_char == 4): #OC03052018
                        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                        srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                        srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0) #OC06052018
                        if((resStokes2 is not None) and (meshRes2 is not None)): #OC30052017
                            #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                            srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                            srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0) #OC06052018
                    else:
                        srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                        if((resStokes2 is not None) and (meshRes2 is not None)): #OC30052017
                            srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))

                    #DEBUG
                    #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)
                    #END DEBUG

                    #MR01112016: write the status of the simulation:  
                    srwl_uti_save_stat_wfr_emit_prop_multi_e(i + 1, total_num_of_particles, filename=log_path)  
                    
                    #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                    
                    #sys.exit(0)
                    iSave = 0

    elif((rank == 0) and (nProc > 1)):

        #nRecv = int(nPartPerProc*nProc/_n_part_avg_proc + 1e-09)
        nRecv = nSentPerProc*(nProc - 1) #Total number of sending acts to be made by all worker processes, and to be received by master

        #print('DEBUG MESSAGE: Actual number of macro-electrons:', nRecv*_n_part_avg_proc)
        
        #DEBUG
        #srwl_uti_save_text("nRecv: " + str(nRecv) + " nPartPerProc: " + str(nPartPerProc) + " nProc: " + str(nProc) + " _n_part_avg_proc: " + str(_n_part_avg_proc), _file_path + ".00.dbg")
        #END DEBUG

        if(resStokes is None):
            resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
        if(workStokes is None):
            workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

        if(_char == 4): #OC30052017 #Cuts of Mutual Intensity vs X & Y
            if(resStokes2 is None):
                resStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
            if(arAuxResSt12 is None):
                lenArResSt12 = len(resStokes.arS) + len(resStokes2.arS)
                arAuxResSt12 = array('f', [0]*lenArResSt12)

            if(workStokes2 is None):
                workStokes2 = SRWLStokes(1, 'f', meshRes2.eStart, meshRes2.eFin, meshRes2.ne, meshRes2.xStart, meshRes2.xFin, meshRes2.nx, meshRes2.yStart, meshRes2.yFin, meshRes2.ny, doMutual)
            if(arAuxWorkSt12 is None):
                lenArWorkSt12 = len(workStokes.arS) + len(workStokes2.arS)
                arAuxWorkSt12 = array('f', [0]*lenArWorkSt12)

        if(_char == 40): #OC03052018 #Intensity and Cuts of Mutual Intensity vs X & Y
            if(resStokes2 is None):
                resStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
            if(resStokes3 is None):
                resStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
            if(arAuxResSt123 is None):
                lenArResSt123 = len(resStokes.arS) + len(resStokes2.arS) + len(resStokes3.arS)
                arAuxResSt123 = array('f', [0]*lenArResSt123)

            if(workStokes2 is None):
                workStokes2 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, _y0, _y0, 1, _mutual=1)
            if(workStokes3 is None):
                workStokes3 = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, _x0, _x0, 1, meshRes.yStart, meshRes.yFin, meshRes.ny, _mutual=1)
            if(arAuxWorkSt123 is None):
                lenArWorkSt123 = len(workStokes.arS) + len(workStokes2.arS) + len(workStokes3.arS)
                arAuxWorkSt123 = array('f', [0]*lenArWorkSt123)

        for i in range(nRecv): #loop over messages from workers

            #DEBUG
            #srwl_uti_save_text("Preparing to receiving # " + str(i), _file_path + ".br.dbg")
            #END DEBUG
           
            #comMPI.Recv([workStokes.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 
            if(workStokes2 is None): #OC31052017
                comMPI.Recv([workStokes.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 
            #else:
            elif(workStokes3 is None): #OC03052018
                comMPI.Recv([arAuxWorkSt12, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 

                lenArWorkSt1 = len(workStokes.arS)
                for i1 in range(lenArWorkSt1): workStokes.arS[i1] = arAuxWorkSt12[i1]
                lenArWorkSt2 = len(workStokes2.arS)
                for i2 in range(lenArWorkSt2): workStokes2.arS[i2] = arAuxWorkSt12[i2 + lenArWorkSt1]

            else: #OC03052018
                comMPI.Recv([arAuxWorkSt123, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 

                lenArWorkSt1 = len(workStokes.arS)
                for i1 in range(lenArWorkSt1): workStokes.arS[i1] = arAuxWorkSt123[i1]
                lenArWorkSt2 = len(workStokes2.arS)
                for i2 in range(lenArWorkSt2): workStokes2.arS[i2] = arAuxWorkSt123[i2 + lenArWorkSt1]
                lenArWorkSt3 = len(workStokes3.arS)
                lenArWorkSt12 = lenArWorkSt1 + lenArWorkSt2
                for i3 in range(lenArWorkSt3): workStokes3.arS[i3] = arAuxWorkSt123[i3 + lenArWorkSt12]

            #if((_char == 4) and (workStokes2 is not None)): #OC30052017 #Cuts of Mutual Intensity vs X & Y
            #    comMPI.Recv([workStokes2.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive 

            #MR20160907 #Save .log and .json files:
            particle_number = (i + 1) * _n_part_avg_proc
            #srwl_uti_save_stat_wfr_emit_prop_multi_e(particle_number, total_num_of_particles)
            srwl_uti_save_stat_wfr_emit_prop_multi_e(particle_number, total_num_of_particles, filename=log_path)

            #DEBUG
            #srwl_uti_save_text("Received intensity # " + str(i), _file_path + ".er.dbg")
            #END DEBUG

            #resStokes.avg_update_same_mesh(workStokes, i + 1)
            #resStokes.avg_update_same_mesh(workStokes, i + 1, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
            multFinAvg = 1 if(_n_part_avg_proc > 1) else ePhIntegMult #OC120714 fixed: the normalization may have been already applied at the previous avaraging on each worker node!

            #print('resStokes.avg_update_same_mesh ... ', end='') #DEBUG
            #t0 = time.time(); #DEBUG

            #DEBUG (test save at i = 0)
            #if(i == 0):
            #    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual)
            #END DEBUG

            #resStokes.avg_update_same_mesh(workStokes, i + 1, 1, multFinAvg) #in the future treat all Stokes components / Polarization, not just s0!
            #resStokes.avg_update_same_mesh(workStokes, i, 1, multFinAvg) #OC15012017 #in the future treat all Stokes components / Polarization, not just s0!
            resStokes.avg_update_same_mesh(workStokes, i, numComp, multFinAvg) #OC15012017 #in the future treat all Stokes components / Polarization, not just s0!

            if((resStokes2 is not None) and (workStokes2 is not None)):
                resStokes2.avg_update_same_mesh(workStokes2, i, numComp, multFinAvg) #OC30052017 #in the future treat all Stokes components / Polarization, not just s0!

            if((resStokes3 is not None) and (workStokes3 is not None)):
                resStokes3.avg_update_same_mesh(workStokes3, i, numComp, multFinAvg) #OC03052018 #in the future treat all Stokes components / Polarization, not just s0!

            #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
            #DEBUG
            #srwl_uti_save_text("Updated Stokes after receiving intensity # " + str(i), _file_path + "." + str(i) + "er.dbg")
            #END DEBUG

            iSave += 1
            if(iSave == _n_save_per):
                #Saving results from time to time in the process of calculation

                #print('srwl_uti_save_intens_ascii ... ', end='') #DEBUG
                #t0 = time.time(); #DEBUG
                
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC26042016
                #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017

                if(_char == 40): #OC03052018
                    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0)
                    if(resStokes2 is not None):
                        #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) 
                        srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC060502018
                        srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                    if(resStokes3 is not None):
                        #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) 
                        srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC060502018
                        srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                elif(_char == 4): #OC03052018
                    #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC060502018
                    srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
                    if((resStokes2 is not None) and (meshRes2 is not None)):
                        #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) 
                        srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC060502018
                        srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0)
                else:
                    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC30052017
                    if((resStokes2 is not None) and (meshRes2 is not None)): #OC30052017
                        srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) 

                #DEBUG
                #srwl_uti_save_intens_ascii(workStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017
                #END DEBUG
                #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
                
                iSave = 0
    #DEBUG
    #srwl_uti_save_text("Exiting srwl_wfr_emit_prop_multi_e", _file_path + "." + str(rank) + "e.dbg")
    #END DEBUG

    if((rank == 0) or (nProc == 1)):
        #Saving final results:
        #if(_file_path != None):
        #if(file_path1 is not None): #OC30052017

        #print('srwl_uti_save_intens_ascii ... ', end='') #DEBUG
        #t0 = time.time(); #DEBUG
            
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC26042016
        #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual) #OC16012017

        if(_char == 40): #OC03052018
            srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 0)
            if(resStokes2 is not None):
                #srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                srwl_uti_save_intens_ascii(resStokes2.arS, resStokes2.mesh, file_path1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), resStokes2.mesh, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
            if(resStokes3 is not None):
                #srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1)
                srwl_uti_save_intens_ascii(resStokes3.arS, resStokes3.mesh, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 1) #OC06052018
                srwl_uti_save_intens_ascii(resStokes3.to_deg_coh(), resStokes3.mesh, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
        elif(_char == 4): #OC03052018
            #srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSaveMutualHorCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
            srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
            srwl_uti_save_intens_ascii(resStokes.to_deg_coh(), meshRes, file_path_deg_coh1, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = 1, _cmplx = 0)
            if((resStokes2 is not None) and (meshRes2 is not None)):
                #srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSaveMutualVerCut, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0))
                srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC06052018
                srwl_uti_save_intens_ascii(resStokes2.to_deg_coh(), meshRes2, file_path_deg_coh2, _n_stokes = 1, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = 0)
        else:
            srwl_uti_save_intens_ascii(resStokes.arS, meshRes, file_path1, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC16012017
            if((resStokes2 is not None) and (meshRes2 is not None)): #OC03052018
                srwl_uti_save_intens_ascii(resStokes2.arS, meshRes2, file_path2, numComp, _arLabels = resLabelsToSave, _arUnits = resUnitsToSave, _mutual = doMutual, _cmplx = (1 if doMutual else 0)) #OC16012017

        #print('completed (lasted', round(time.time() - t0, 6), 's)') #DEBUG
        #return resStokes
        if(resStokes2 is None): #OC30052017
            return resStokes
        #else:
        elif(resStokes3 is None): #OC03052018
            return resStokes, resStokes2
        else:
            return resStokes, resStokes2, resStokes3
    else:
        return None

#****************************************************************************
#****************************************************************************
#Import of modules requiring classes defined in this smodule
#****************************************************************************
#****************************************************************************
from srwl_uti_src import *

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
       _inSig[0]: RMS size of teh Gaussian in first dimension
       _inSig[1]: (optional) RMS size of a 2D Gaussian in second dimension
       _inSig[2]: (optional) coefficient before cross-term in exponent argument of a 2D Gaussian
       i.e. _inSig[] = [sigX, sigY, alp} defines a "tilted" normalized 2D Gaussian (vs x, y): 
       (sqrt(1 - (alp*sigX*sigY)**2)/(2*Pi*sigX*sigY))*exp(-x**2/(2*sigX**2) - y**2/(2*sigY^2) - alp*x*y)
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
