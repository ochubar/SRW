#############################################################################
# SRWLib for Python v 0.11
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
import srwlpy as srwl
from array import *
from math import *
from copy import *
import random
import sys
import os
import traceback
import uti_math
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

class SRWLMagFldM(SRWLMagFld):
    """Magnetic Field: Multipole Magnet"""
    
    def __init__(self, _G=0, _m=2, _n_or_s='n', _Leff=0, _Ledge=0):
        """
        :param _G: field parameter [T] for dipole, [T/m] for quadrupole (negative means defocusing for x), [T/m^2] for sextupole, [T/m^3] for octupole
        :param _m: multipole order 1 for dipole, 2 for quadrupoole, 3 for sextupole, 4 for octupole
        :param _n_or_s: normal ('n') or skew ('s')
        :param _Leff: effective length [m]
        :param _Ledge: "soft" edge length for field variation from 10% to 90% [m]; G/(1 + ((z-zc)/d)^2)^2 fringe field dependence is assumed
        """
        self.G = _G
        self.m = _m
        self.n_or_s = _n_or_s
        self.Leff = _Leff
        self.Ledge = _Ledge

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
        self.arHarm = [SRWLMagFldH()]*_nHarm

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
    
    def __init__(self, _arMagFld=None, _arXc=None, _arYc=None, _arZc=None):
        """
        :param _arMagFld: magnetic field structures array
        :param _arXc: horizontal center positions of magnetic field elements in arMagFld array [m]
        :param _arYc: vertical center positions of magnetic field elements in arMagFld array [m]
        :param _arZc: longitudinal center positions of magnetic field elements in arMagFld array [m]
        """
        self.arMagFld = [] if _arMagFld is None else _arMagFld
        self.arXc = array('d') if _arXc is None else _arXc
        self.arYc = array('d') if _arYc is None else _arYc
        self.arZc = array('d') if _arZc is None else _arZc

    def allocate(self, _nElem):
        self.arMagFld = [SRWLMagFld()]*_nElem
        self.arXc = array('d', [0]*_nElem)
        self.arYc = array('d', [0]*_nElem)
        self.arZc = array('d', [0]*_nElem)

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
        self.arX = array('d') if _arX is None else _arX
        self.arY = array('d') if _arY is None else _arY
        self.arZ = array('d') if _arZ is None else _arZ
        self.arXp = array('d') if _arXp is None else _arXp
        self.arYp = array('d') if _arYp is None else _arYp
        self.arZp = array('d') if _arZp is None else _arZp
        if _arBx != None: #by default, arBx, _arBy, arBz is not created
            self.arBx = _arBx
        if _arBy != None:
            self.arBy = _arBy
        if _arBz != None:
            self.arBz = _arBz
        self.np = _np
        self.ctStart = _ctStart
        self.ctEnd = _ctEnd
        self.partInitCond = SRWLParticle() if _partInitCond is None else _partInitCond

    def allocate(self, _np, _allB=False):
        _np = int(_np)
        self.arX = array('d', [0] * _np)
        self.arXp = array('d', [0] * _np)
        self.arY = array('d', [0] * _np)
        self.arYp = array('d', [0] * _np)
        self.arZ = array('d', [0] * _np)
        self.arZp = array('d', [0] * _np)
        self.np = _np
        if _allB == True:
            self.arBx = array('d', [0] * _np)
            self.arBy = array('d', [0] * _np)
            self.arBz = array('d', [0] * _np)

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
    """Gaussian Beam"""
    
    def __init__(self, _x=0, _y=0, _z=0, _xp=0, _yp=0, _avgPhotEn=1, _pulseEn=1, _repRate=1, _polar=1, _sigX=10e-06,
                 _sigY=10e-06, _sigT=1e-15, _mx=0, _my=0):
        """
        :param _x: average horizontal coordinates of waist [m]
        :param _y: average vertical coordinates of waist [m]
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
        if(_mesh.arSurf != None):
            try:
                lenArSurf = len(_mesh.arSurf)
                if(lenArSurf > 0):
                    self.arSurf = array('d', [0]*lenArSurf)
                    for i in range(lenArSurf):
                        self.arSurf[i] = _mesh.arSurf[i]
            except:
                pass

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
            nTot *= nTot
        nTot *= 4 #array length of all Stokes components
        #eventually allow for storage of less than 4 Stokes components!
        
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

        nTot = _n_comp*self.mesh.ne*self.mesh.nx*self.mesh.ny #eventually allow for storage of less than 4 Stokes components!
        if(self.mutual > 0): nTot *= nTot

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
            nStPt *= nStPt
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

        perE = eNpWfr
        perXp = perE*eNpWfr
        perX = perXp*xNpWfr
        perYp = perX*xNpWfr
        perY = perYp*yNpWfr
        nRadWfr = perY*yNpWfr
        
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

                                    fInterp = 0
                                    if(not(doZeroFy or doZeroFyp or doZeroFx or doZeroFxp or doZeroFe or doZeroFep)):
                                        a000000 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f100000 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f010000 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f001000 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f000100 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        f000010 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        f000001 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]

                                        f110000 = _more_stokes.arS[iOfstSt + iep1 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f101000 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f100100 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        f100010 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        f100001 = _more_stokes.arS[iOfstSt + iep1 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]

                                        f011000 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy0_perY]
                                        f010100 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        f010010 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        f010001 = _more_stokes.arS[iOfstSt + iep0 + ie1_perE + ixp0_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]
                                   
                                        f001100 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix1_perX_p_iyp0_perYp_p_iy0_perY]
                                        f001010 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp1_perYp_p_iy0_perY]
                                        f001001 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp1_perXp_p_ix0_perX_p_iyp0_perYp_p_iy1_perY]

                                        f000110 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp1_perYp_p_iy0_perY]
                                        f000101 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix1_perX_p_iyp0_perYp_p_iy1_perY]

                                        f000011 = _more_stokes.arS[iOfstSt + iep0 + ie0_perE + ixp0_perXp_p_ix0_perX_p_iyp1_perYp_p_iy1_perY]

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

                                    #self.arS[ir] = (self.arS[ir]*_iter + fInterp)/(_iter + 1)
                                    self.arS[ir] = (self.arS[ir]*_iter + _mult*fInterp)/(_iter + 1)
                                    
                                    ir += 1
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
    #unitElFld = 1 #electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2) ?
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

    def allocate(self, _ne, _nx, _ny, EXNeeded=1, EYNeeded=1, typeE='f'):
        #print('') #debugging
        #print('          (re-)allocating: old point numbers: ne=',self.mesh.ne,' nx=',self.mesh.nx,' ny=',self.mesh.ny) #,' type:',self.numTypeElFld)
        #print('                           new point numbers: ne=',_ne,' nx=',_nx,' ny=',_ny) #,' type:',typeE)
        nTot = 2*_ne*_nx*_ny #array length to store one component of complex electric field
        nMom = 11*_ne
        if EXNeeded:
            #print('          trying to (re-)allocate Ex ... ', end='')  
            del self.arEx
            self.arEx = array(typeE, [0]*nTot)
            #print('done')           
            if len(self.arMomX) != nMom:
                del self.arMomX
                self.arMomX = array('d', [0]*nMom)
        if EYNeeded:
            #print('          trying to (re-)allocate Ey ... ', end='')  
            #del self.arEy
            self.arEy = array(typeE, [0]*nTot)
            #print('done')
            if len(self.arMomY) != nMom:
                del self.arMomY
                self.arMomY = array('d', [0]*nMom)
        self.numTypeElFld = typeE
        self.mesh.ne = _ne
        self.mesh.nx = _nx
        self.mesh.ny = _ny

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

    def calc_stokes(self, _stokes):
        """Calculate Stokes parameters from Electric Field"""
        if(_stokes.mutual <= 0):
            nTot = self.mesh.ne*self.mesh.nx*self.mesh.ny
            #if(type(_stokes).__name__ != 'SRWLStokes')):
            if(isinstance(_stokes, SRWLStokes) == False):
                raise Exception("Incorrect Stokes parameters object submitted") 
            nTotSt = nTot*4
            nTot2 = nTot*2
            nTot3 = nTot*3
            if(_stokes.arS != None):
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
                _stokes.arS[i + nTot] = intLinX - intLinY #check sign
                _stokes.arS[i + nTot2] = -2*(reEx*reEy + imEx*imEy) #check sign
                _stokes.arS[i + nTot3] = 2*(-reEx*reEy + imEx*imEy) #check sign
            _stokes.mesh.set_from_other(self.mesh)
            
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

            nTot = eNpRes*xNpRes*yNpRes
            nTot1 = nTot*nTot
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
            auxArEx = array('f', [0]*nTotAux)
            auxArEy = array('f', [0]*nTotAux)
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

                                    intLinX = reEx*reExT + imEx*imExT
                                    intLinY = reEy*reEyT + imEy*imEyT#; print(intLinX, intLinY)
                                    _stokes.arS[ir] = intLinX + intLinY
                                    _stokes.arS[ir + nTot1] = intLinX - intLinY #check sign
                                    _stokes.arS[ir + nTot2] = -reEx*reEyT - reExT*reEy - imEx*imEyT - imExT*imEy #-2*(reEx*reEy + imEx*imEy) #check sign
                                    _stokes.arS[ir + nTot3] = -reEx*reEyT - reExT*reEy + imEx*imEyT + imExT*imEy #2*(-reEx*reEy + imEx*imEy) #check sign
                                    ir += 1

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
        :param _nx: number of transmission data points in the horizontaldirection
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
        if((_arTr == None) or ((len(_arTr) != _ne*_nx*_ny*2) and (_ne*_nx*_ny > 0))):
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

class SRWLOptMirPl(SRWLOptMir):
    """Optical Element: Mirror: Plane"""
    
    def __init__(self, 
                 _size_tang=1, _size_sag=1, _ap_shape='r', _sim_meth=2, _npt=100, _nps=100, _treat_in_out=1, _ext_in=0, _ext_out=0,
                 _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0, _x=0, _y=0,
                 _refl=1, _n_ph_en=1, _n_ang=1, _n_comp=1, _ph_en_start=1000., _ph_en_fin=1000., _ph_en_scale_type='lin', _ang_start=0, _ang_fin=0, _ang_scale_type='lin'):
        """
        :param _r_sag: sagital radius of curvature at mirror center [m]
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
        self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

class SRWLOptMirEl(SRWLOptMir):
    """Optical Element: Mirror: Ellipsoid"""
    
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
        self.set_dim_sim_meth(_size_tang, _size_sag, _ap_shape, _sim_meth, _npt, _nps, _treat_in_out, _ext_in, _ext_out)
        self.set_orient(_nvx, _nvy, _nvz, _tvx, _tvy, _x, _y)
        self.set_reflect(_refl, _n_ph_en, _n_ang, _n_comp, _ph_en_start, _ph_en_fin, _ph_en_scale_type, _ang_start, _ang_fin, _ang_scale_type)

#class SRWLOptMirTor(SRWLOptMir):
#    """Optical Element: Mirror: Toroid (to be developed)"""
#    
#    def __init__(self, _rt=1, _rs=1):
#        """
#        :param _rt: tangential (major) radius [m]
#        :param _rs: sagittal (minor) radius [m]
#        """
#        self.radTan = _rt
#        self.radSag = _rs
#        
#        #finishing of the mirror setup requires calling these 3 functions (with their required arguments):
#        self.set_reflect()
#        self.set_dim_sim_meth()
#        self.set_orient()

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
    def __init__(self, _d_sp, _psi0r, _psi0i, _psi_hr, _psi_hi, _psi_hbr, _psi_hbi, _tc, _ang_as, _nvx=0, _nvy=0, _nvz=-1, _tvx=1, _tvy=0):
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
        """
        #"""
        #The Miller Inices are removed from this input (after discussion with A. Suvorov), because _d_sp already incorporates this information:
        #:param _h1: 1st index of diffraction vector (John's hMilND)
        #:param _h2: 2nd index of diffraction vector (John's kMilND)
        #:param _h3: 3rd index of diffraction vector (John's lMilND)
        #However, a member-function may be added here to calculate _d_sp from teh Miller Indices and material constant(s)
        
        #Moved to Propagation Parameters
        #:param _sx: horizontal coordinate of optical axis after crystal [m] (John's )
        #:param _sy: vertical coordinate of optical axis after crystal [m] (John's )
        #:param _sz: longitudinal coordinate of optical axis after crystal [m] (John's )

        #to go to member functions (convenience derived parameters)
        #:param _aChi: crystal roll angle (John's )
        #:param _aPsi: crystal yaw angle (John's )
        #:param _aThe: crystal theta angle (John's )
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
        """
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

        tolAng = 1.e-06
        if(abs(_ang_dif_pl) < tolAng): #case of the vertical deflection plane
            return [[tv, sv, nv], [rx, ry, rz]]
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

            ex = None; ey = None
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
            return [[prodMV(mr, tv), prodMV(mr, sv), prodMV(mr, nv)], [ex, ey, ez]]

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
            [6]: Horizontal Resolution modification factor at Resizing
            [7]: Vertical Range modification factor at Resizing
            [8]: Vertical Resolution modification factor at Resizing
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
        if(_arOpt == None):
            self.arOpt = []
        self.arProp = _arProp #list of lists of propagation parameters to be used for individual optical elements
        if(_arProp == None):
            self.arProp = []
            
    def allocate(self, _nElem):
        self.arOpt = [SRWLOpt()]*_nElem
        self.arProp = [[0]*17]*_nElem  

#****************************************************************************
#****************************************************************************
#Setup some transmission-type optical elements
#****************************************************************************
#****************************************************************************
def srwl_opt_setup_CRL(_foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n, _wall_thick, _xc, _yc, _void_cen_rad=None, _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
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
    :return: transmission (SRWLOptT) type optical element which simulates CRL
    """
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
                    return sectLen - 2*sqrt(radE2 - rE2)
            elif(rE2 < radE2):
                return sectLen - 2*sqrt(radE2 - rE2)
        return sectLen

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

    foc_len = (0.5*_r_min/(_n*arDelta[int(0.5*ne)]))
    print('Optical Element Setup: CRL Focal Length:', foc_len, 'm')

    fx = 1e+23
    fy = 1e+23
    if(_foc_plane != 1):
        fy = foc_len
    if(_foc_plane != 2):
        fx = foc_len
    
    opT = SRWLOptT(nx, ny, rx, ry, None, 1, fx, fy, _xc, _yc, ne, _e_start, _e_fin)

    #print(ne, _e_start, _e_fin)

    halfApert = 0.5*_apert_h
    halfApertV = 0.5*_apert_v

    if(_foc_plane == 2): #1D lens, vertical is focusing plane
        halfApert = halfApertV
    elif(_foc_plane == 3): #2D lens
        if(halfApert > halfApertV):
            halfApert = halfApertV
    
    hx = rx/(nx - 1)
    hy = ry/(ny - 1)

    #Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    ofst = 0
    y = -0.5*ry #CRL is always centered on the grid, however grid can be shifted
    for iy in range(ny):
        x = -0.5*rx
        for ix in range(nx):
            pathInBody = _n*ray_path_in_one_CRL(x, y, _foc_plane, _shape, halfApert, _r_min, _wall_thick)
    
            if(_void_cen_rad != None): #eventually subtract path in voids
                pathInBody -= ray_path_in_spheres(x, y, _void_cen_rad)

            for ie in range(ne):
                opT.arTr[ofst] = exp(-0.5*pathInBody/arAttenLen[ie]) #amplitude transmission
                opT.arTr[ofst + 1] = -arDelta[ie]*pathInBody #optical path difference
                ofst += 2
            x += hx
        y += hy
        
    return opT

#****************************************************************************
def srwl_opt_setup_cyl_fiber(_foc_plane, _delta_ext, _delta_core, _atten_len_ext, _atten_len_core, _diam_ext, _diam_core, _xc, _yc):
    """
    Setup Transmission type Optical Element which simulates Cylindrical Fiber
    :param _foc_plane: plane of focusing: 1- horizontal (i.e. fiber is parallel to vertical axis), 2- vertical (i.e. fiber is parallel to horizontal axis)
    :param _delta_ext: refractive index decrement of extenal layer
    :param _delta_core: refractive index decrement of core
    :param _atten_len_ext: attenuation length [m] of external layer
    :param _atten_len_core: attenuation length [m] of core
    :param _diam_ext: diameter [m] of external layer
    :param _diam_core: diameter [m] of core
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :return: transmission (SRWLOptT) type optical element which simulates Cylindrical Fiber
    """

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
    return opT

#****************************************************************************
def srwl_opt_setup_surf_height_1d(_height_prof_data, _dim, _ang, _ang_r=0, _amp_coef=1, _ar_arg_long=None, _nx=0, _ny=0, _size_x=0, _size_y=0):
    """
    Setup Transmission type optical element with 1D (mirror or grating) surface Heght Profile data
    :param _height_prof_data: two- or one-column table containing, in case of two columns: longitudinal position in [m] (1st column) and the Height Profile in [m] (2nd column) data; in case of one column, it contains the Height Profile data
    :param _dim: orientation of the reflection (deflection) plane; can be 'x' or 'y'
    :param _ang: grazing angle (between input optical axis and mirror/grating plane)
    :param _ang_r: reflection angle (between output optical axis and mirror/grating plane)
    :param _amp_coef: height profile "amplification coefficient"
    :param _ar_arg_long: optional array of longitudinal position (along mirror/grating) in [m]; if _ar_arg_long != None, any longitudinal position contained in _height_prof_data is ignored
    :param _nx: optional number of points in horizontal dimension of the output transmission optical element
    :param _ny: optional number of points in vertical dimension of the output transmission optical element
    :param _size_x: optional horizontal transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :param _size_y: optional vertical transverse size of the output transmission optical element (if <=0: _height_prof_data, _dim, _ar_arg_long, _ar_arg_tr data is used)
    :return: transmission (SRWLOptT) type optical element which simulates the effect of surface height error
    """
    #To test all options!

    if(_ang_r == 0): _ang_r = _ang
    sinAng = sin(_ang)
    sinAngR = sin(_ang_r)

    if _ar_arg_long == None:
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

    optSlopeErr = SRWLOptT(nx, ny, sizeX, sizeY)

    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)

    y = auxMesh.yStart
    hApprox = 0
    ipStart = 0
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

        #x = optSlopeErr.x - 0.5*optSlopeErr.rx
        x = auxMesh.xStart
        
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

            optSlopeErr.arTr[ofst] = 1. #Amplitude Transmission
            optSlopeErr.arTr[ofst + 1] = 0. #Optical Path Difference
            if(hApprox != 0):
                #optSlopeErr.arTr[ofst + 1] = -2*sinAng*hApprox #Optical Path Difference (to check sign!)
                optSlopeErr.arTr[ofst + 1] = -(sinAng + sinAngR)*hApprox*_amp_coef #Optical Path Difference (to check sign!)
                #print(ix, iy, optSlopeErr.arTr[ofst + 1])
            x += xStep
        y += yStep
    return optSlopeErr

#****************************************************************************
def srwl_opt_setup_surf_height_2d(_height_prof_data, _dim, _ang, _ang_r=0, _amp_coef=1, _ar_arg_long=None, _ar_arg_tr=None, _nx=0, _ny=0, _size_x=0, _size_y=0):
    """
    Setup Transmission type optical element with 2D (mirror or grating) surface Heght Profile data
    :param _height_prof_data: a matrix (2D array) containing the Height Profile data in [m]; if _ar_height_prof_x==None and _ar_height_prof_y==None: the first column in _height_prof_data is assumed to be the "longitudinal" position [m] and first row the "transverse" position [m], and _height_prof_data[0][0] is not used; otherwise the "longitudinal" and "transverse" positions on the surface are assumed to be given by _ar_height_prof_x, _ar_height_prof_y 
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

    if(_ang_r == 0): _ang_r = _ang
    sinAng = sin(_ang)
    sinAngR = sin(_ang_r)
    
    #argHeightProfData = _ar_arg_long
    if _ar_arg_long == None:
        npData = len(_height_prof_data[0]) - 1
        sizeLong = _height_prof_data[0][npData - 1] - _height_prof_data[0][1]
    else:
        npData = len(_ar_arg_long)
        sizeLong = _ar_arg_long[npData - 1] - _ar_arg_long[0]
        
    sizeLongProj = sizeLong*sinAngR

    if _ar_arg_tr == None:
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

    print(auxMesh.xStart, auxMesh.xFin, auxMesh.nx, xStep)
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
            if _ar_arg_long == None: y1 = _height_prof_data[0][ipStart]*sinAngR
            else: y1 = _ar_arg_long[ipStart - 1]*sinAngR
            
            for i in range(ipStart + 1, npData + 1):
            #for i in range(ipStart + 1, npData):
                #y2 = argHeightProfData[i]*sinAngR
                if _ar_arg_long == None: y2 = _height_prof_data[0][i]*sinAngR
                else: y2 = _ar_arg_long[i - 1]*sinAngR
                
                if((y1 <= y) and (y < y2)):
                    ipStart = i - 1
                    break
                y1 = y2

        elif('x' in _dim):
            ipStart = 1
            if _ar_arg_tr == None: y1 = _height_prof_data[ipStartTr][0]
            else: y1 = _ar_arg_tr[ipStartTr - 1]
            
            for i in range(ipStartTr + 1, npDataTr + 1):
            #for i in range(ipStartTr + 1, npDataTr):
                if _ar_arg_tr == None: y2 = _height_prof_data[i][0]
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
              
                if _ar_arg_tr == None: x1 = _height_prof_data[ipStartTr][0]
                else: x1 = _ar_arg_tr[ipStartTr - 1]

                #print(ipStartTr + 1, npDataTr + 1)
                
                for i in range(ipStartTr + 1, npDataTr + 1):
                #for i in range(ipStartTr + 1, npDataTr):
                    if _ar_arg_tr == None: x2 = _height_prof_data[i][0]
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
                if _ar_arg_long == None: x1 = _height_prof_data[0][ipStart]*sinAngR
                else: x1 = _ar_arg_long[ipStart - 1]*sinAngR
                
                for i in range(ipStart + 1, npData + 1):
                #for i in range(ipStart + 1, npData):
                    #x2 = argHeightProfData[i]*sinAngR
                    if _ar_arg_long == None: x2 = _height_prof_data[0][i]*sinAngR
                    else: x2 = _ar_arg_long[i - 1]*sinAngR
                    
                    if((x1 <= x) and (x < x2)):
                        ipStart = i - 1
                        break
                    x1 = x2

            if _ar_arg_long != None: ipStart -= 1
            if _ar_arg_tr != None: ipStartTr -= 1

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
    return optSlopeErr

#****************************************************************************
#****************************************************************************
#Auxiliary utility functions
#****************************************************************************
#****************************************************************************
#Moved to uti_math.py
##def srwl_uti_interp_1d(_x, _x_min, _x_step, _nx, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
##    """
##    Interpolate 1D function value tabulated on equidistant mesh, using polynomial interpolation
##    :param _x: argument at which function value should be calculated
##    :param _x_min: minimal argument value of the tabulated function
##    :param _x_step: step of mesh at which function is tabulated
##    :param _nx: number of points in mesh at which function is tabulated
##    :param _ar_f: tabulated function list or array
##    :param _ord: order of polynomial interpolation (1- linear, 2- quadratic, 3- cubic)
##    :param _ix_per: argument index period of function data alignment (e.g. to interpolate one component of complex data, or in one dimension of multi-dimensional data)
##    :param _ix_ofst: argument index offset of function data alignment
##    :return: function value found by polynomial interpolation
##    """
##    if(_ord == 1):
##        i0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
##        if(i0 < 0):
##            i0 = 0
##        elif(i0 >= _nx - 1):
##            i0 = _nx - 2
##        i1 = i0 + 1
##        f0 = _ar_f[i0*_ix_per + _ix_ofst]
##        f1 = _ar_f[i1*_ix_per + _ix_ofst]
##        t = (_x - (_x_min + _x_step*i0))/_x_step
##        return f0 + (f1 - f0)*t
##    elif(_ord == 2):
##        i0 = int(round((_x - _x_min)/_x_step))
##        if(i0 < 1):
##            i0 = 1
##        elif(i0 >= _nx - 1):
##            i0 = _nx - 2
##        im1 = i0 - 1
##        i1 = i0 + 1
##        t = (_x - (_x_min + _x_step*i0))/_x_step
##        a0 = _ar_f[i0*_ix_per + _ix_ofst]
##        fm1 = _ar_f[im1*_ix_per + _ix_ofst]
##        f1 = _ar_f[i1*_ix_per + _ix_ofst]
##        a1 = 0.5*(f1 - fm1)
##        a2 = 0.5*(fm1 + f1 - 2*a0)
##        return a0 + t*(a1 + t*a2)
##    elif(_ord == 3):
##        i0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
##        if(i0 < 1):
##            i0 = 1
##        elif(i0 >= _nx - 2):
##            i0 = _nx - 3
##        im1 = i0 - 1
##        i1 = i0 + 1
##        i2 = i0 + 2
##        t = (_x - (_x_min + _x_step*i0))/_x_step
##        a0 = _ar_f[i0*_ix_per + _ix_ofst]
##        fm1 = _ar_f[im1*_ix_per + _ix_ofst]
##        f1 = _ar_f[i1*_ix_per + _ix_ofst]
##        f2 = _ar_f[i2*_ix_per + _ix_ofst]
##        a1 = -0.5*a0 + f1 - f2/6. - fm1/3.
##        a2 = -a0 + 0.5*(f1 + fm1)
##        a3 = 0.5*(a0 - f1) + (f2 - fm1)/6.
##        return a0 + t*(a1 + t*(a2 + t*a3))
##    return 0

#****************************************************************************
##def srwl_uti_interp_2d(_x, _y, _x_min, _x_step, _nx, _y_min, _y_step, _ny, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
#Moved to uti_math.py
##    """
##    Interpolate 2D function value tabulated on equidistant rectangular mesh and represented by C-aligned flat array, using polynomial interpolation
##    :param _x: first argument at which function value should be calculated
##    :param _y: second argument at which function value should be calculated
##    :param _x_min: minimal value of the first argument of the tabulated function
##    :param _x_step: step of the first argument at which function is tabulated
##    :param _nx: number of points vs first argument at which function is tabulated
##    :param _y_min: minimal value of the second argument of the tabulated function
##    :param _y_step: step of the second argument at which function is tabulated
##    :param _ny: number of points vs second argument at which function is tabulated
##    :param _ar_f: function tabulated on 2D mesh, aligned as "flat" C-type list or array (first argument is changing most frequently)
##    :param _ord: "order" of polynomial interpolation (1- bi-linear (on 4 points), 2- "bi-quadratic" (on 6 points), 3- "bi-cubic" (on 12 points))
##    :param _ix_per: period of first argument index of the function data alignment (e.g. to interpolate one component of complex data, or in one dimension of multi-dimensional data)
##    :param _ix_ofst: offset of the first argument index in function data alignment
##    :return: function value found by 2D polynomial interpolation
##    """
##    if(_ord == 1): #bi-linear interpolation based on 4 points
##        ix0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
##        if(ix0 < 0):
##            ix0 = 0
##        elif(ix0 >= _nx - 1):
##            ix0 = _nx - 2
##        ix1 = ix0 + 1
##        tx = (_x - (_x_min + _x_step*ix0))/_x_step
##        
##        iy0 = int(trunc((_y - _y_min)/_y_step + 1.e-09))
##        if(iy0 < 0):
##            iy0 = 0
##        elif(iy0 >= _ny - 1):
##            iy0 = _ny - 2
##        iy1 = iy0 + 1
##        ty = (_y - (_y_min + _y_step*iy0))/_y_step
##
##        nx_ix_per = _nx*_ix_per
##        iy0_nx_ix_per = iy0*nx_ix_per
##        iy1_nx_ix_per = iy1*nx_ix_per
##        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
##        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
##        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        a10 = f10 - a00
##        a01 = f01 - a00
##        a11 = a00 - f01 - f10 + f11
##        return a00 + tx*(a10 + ty*a11) + ty*a01
##
##    elif(_ord == 2): #bi-quadratic interpolation based on 6 points
##        ix0 = int(round((_x - _x_min)/_x_step))
##        if(ix0 < 1):
##            ix0 = 1
##        elif(ix0 >= _nx - 1):
##            ix0 = _nx - 2
##        ixm1 = ix0 - 1
##        ix1 = ix0 + 1
##        tx = (_x - (_x_min + _x_step*ix0))/_x_step
##
##        iy0 = int(round((_y - _y_min)/_y_step))
##        if(iy0 < 1):
##            iy0 = 1
##        elif(iy0 >= _ny - 1):
##            iy0 = _ny - 2
##        iym1 = iy0 - 1
##        iy1 = iy0 + 1
##        ty = (_y - (_y_min + _y_step*iy0))/_y_step
##
##        nx_ix_per = _nx*_ix_per
##        iym1_nx_ix_per = iym1*nx_ix_per
##        iy0_nx_ix_per = iy0*nx_ix_per
##        iy1_nx_ix_per = iy1*nx_ix_per
##        ixm1_ix_per_p_ix_ofst = ixm1*_ix_per + _ix_ofst
##        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
##        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
##        fm10 = _ar_f[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
##        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        f0m1 = _ar_f[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        a10 = 0.5*(f10 - fm10)
##        a01 = 0.5*(f01 - f0m1)
##        a11 = a00 - f01 - f10 + f11
##        a20 = 0.5*(f10 + fm10) - a00
##        a02 = 0.5*(f01 + f0m1) - a00
##        return a00 + tx*(a10 + tx*a20 + ty*a11) + ty*(a01 + ty*a02)
##    
##    elif(_ord == 3): #bi-cubic interpolation based on 12 points
##        ix0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
##        if(ix0 < 1):
##            ix0 = 1
##        elif(ix0 >= _nx - 2):
##            ix0 = _nx - 3
##        ixm1 = ix0 - 1
##        ix1 = ix0 + 1
##        ix2 = ix0 + 2
##        tx = (_x - (_x_min + _x_step*ix0))/_x_step
##
##        iy0 = int(trunc((_y - _y_min)/_y_step + 1.e-09))
##        if(iy0 < 1):
##            iy0 = 1
##        elif(iy0 >= _ny - 2):
##            iy0 = _ny - 3
##        iym1 = iy0 - 1
##        iy1 = iy0 + 1
##        iy2 = iy0 + 2
##        ty = (_y - (_y_min + _y_step*iy0))/_y_step
##
##        nx_ix_per = _nx*_ix_per
##        iym1_nx_ix_per = iym1*nx_ix_per
##        iy0_nx_ix_per = iy0*nx_ix_per
##        iy1_nx_ix_per = iy1*nx_ix_per
##        iy2_nx_ix_per = iy2*nx_ix_per
##        ixm1_ix_per_p_ix_ofst = ixm1*_ix_per + _ix_ofst
##        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
##        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
##        ix2_ix_per_p_ix_ofst = ix2*_ix_per + _ix_ofst
##        f0m1 = _ar_f[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f1m1 = _ar_f[iym1_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        fm10 = _ar_f[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
##        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        f20 = _ar_f[iy0_nx_ix_per + ix2_ix_per_p_ix_ofst]
##        fm11 = _ar_f[iy1_nx_ix_per + ixm1_ix_per_p_ix_ofst]
##        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        f21 = _ar_f[iy1_nx_ix_per + ix2_ix_per_p_ix_ofst]
##        f02 = _ar_f[iy2_nx_ix_per + ix0_ix_per_p_ix_ofst]
##        f12 = _ar_f[iy2_nx_ix_per + ix1_ix_per_p_ix_ofst]
##        a10 = -0.5*a00 + f10 - f20/6 - fm10/3
##        a01 = -0.5*a00 + f01 - f02/6 - f0m1/3
##        a11 = -0.5*(f01 + f10) + (f02 - f12 + f20 - f21)/6 + (f0m1 - f1m1 + fm10 - fm11)/3 + f11
##        a20 = -a00 + 0.5*(f10 + fm10)
##        a02 = -a00 + 0.5*(f01 + f0m1)
##        a21 = a00 - f01 + 0.5*(f11 - f10 - fm10 + fm11)
##        a12 = a00 - f10 + 0.5*(f11 - f01 - f0m1 + f1m1)
##        a30 = 0.5*(a00 - f10) + (f20 - fm10)/6
##        a03 = 0.5*(a00 - f01) + (f02 - f0m1)/6
##        a31 = 0.5*(f01 + f10 - f11 - a00) + (f21 + fm10 - f20 - fm11)/6
##        a13 = 0.5*(f10 - f11 - a00 + f01) + (f0m1 + f12 - f02 - f1m1)/6
##        return a00 + tx*(a10 + tx*(a20 + tx*(a30 + ty*a31) + ty*a21) + ty*a11) + ty*(a01 + ty*(a02 + ty*(a03 + tx*a13) + tx*a12))
##    return 0

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
    elif _out_u == 'um': x = (_Light_eV_mu*1.e-03)/x #this had to be modofoed because of non-ascii "mu" symbol that did not compile on Py 2.7
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
        #resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], -1) #MPI module dynamic load
        resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], 0) #MPI module dynamic load
        #multiply re-import won't hurt; but it would be better to avoid this(?)
        MPI = resImpMPI4Py.MPI
        comMPI = MPI.COMM_WORLD
        rankMPI = comMPI.Get_rank()
        if(rankMPI == 0):
            return True
        else:
            return False
    except:
        return True

#**********************Auxiliary function to write tabulated resulting Intensity data to an ASCII file:
def srwl_uti_save_intens_ascii(_ar_intens, _mesh, _file_path, _n_stokes=1, _arLabels=['Photon Energy', 'Horizontal Position', 'Vertical Position', 'Intensity'], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'], _mutual=0):
    f = open(_file_path, 'w')
    arLabelUnit = [_arLabels[i] + ' [' + _arUnits[i] + ']' for i in range(4)]
    f.write('#' + arLabelUnit[3] + ' (C-aligned, inner loop is vs ' + _arLabels[0] + ', outer loop vs ' + _arLabels[2] + ')\n')
    f.write('#' + repr(_mesh.eStart) + ' #Initial ' + arLabelUnit[0] + '\n')
    f.write('#' + repr(_mesh.eFin) + ' #Final ' + arLabelUnit[0] + '\n')
    f.write('#' + repr(_mesh.ne) + ' #Number of points vs ' + _arLabels[0] + '\n')
    f.write('#' + repr(_mesh.xStart) + ' #Initial ' + arLabelUnit[1] + '\n')
    f.write('#' + repr(_mesh.xFin) + ' #Final ' + arLabelUnit[1] + '\n')
    f.write('#' + repr(_mesh.nx) + ' #Number of points vs ' + _arLabels[1] + '\n')
    f.write('#' + repr(_mesh.yStart) + ' #Initial ' + arLabelUnit[2] + '\n')
    f.write('#' + repr(_mesh.yFin) + ' #Final ' + arLabelUnit[2] + '\n')
    f.write('#' + repr(_mesh.ny) + ' #Number of points vs ' + _arLabels[2] + '\n')
    nComp = 1
    if _n_stokes > 0:
        f.write('#' + repr(_n_stokes) + ' #Number of components\n')
        nComp = _n_stokes
    nRadPt = _mesh.ne*_mesh.nx*_mesh.ny
    if(_mutual > 0): nRadPt *= nRadPt
    
    nVal = nRadPt*nComp #_mesh.ne*_mesh.nx*_mesh.ny*nComp
    for i in range(nVal): #write all data into one column using "C-alignment" as a "flat" 1D array
        f.write(' ' + repr(_ar_intens[i]) + '\n')
    f.close()

#**********************Auxiliary function to write auxiliary/debugging information to an ASCII file:
def srwl_uti_save_text(_text, _file_path):
    f = open(_file_path, 'w')
    f.write(_text + '\n')
    f.close()

#**********************Auxiliary function to read-in data comumns from ASCII file (2D table):
def srwl_uti_read_data_cols(_file_path, _str_sep, _i_col_start=0, _i_col_end=-1, _n_line_skip=0):
    """
    Auxiliary function to read-in data comumns from ASCII file (2D table)
    :param _file_path: full path (including file name) to the file
    :param _i_col_start: initial data column to read
    :param _i_col_end: final data column to read
    :param _str_sep: column separation symbol(s) (string)
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
def srwl_wfr_emit_prop_multi_e(_e_beam, _mag, _mesh, _sr_meth, _sr_rel_prec, _n_part_tot, _n_part_avg_proc=1, _n_save_per=100, _file_path=None, _sr_samp_fact=-1, _opt_bl=None, _pres_ang=0, _char=0, _x0=0, _y0=0, _e_ph_integ=0, _rand_meth=1):
    """
    Calculate Stokes Parameters of Emitted (and Propagated, if beamline is defined) Partially-Coherent SR
    :param _e_beam: Finite-Emittance e-beam (SRWLPartBeam type)
    :param _mag: Magnetic Field container (magFldCnt type)
    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
    :param _sr_meth: SR Electric Field calculation method to be used (0- "manual", 1- "auto-undulator", 2- "auto-wiggler")
    :param _sr_rel_prec: relative precision for SR Electric Field calculation (usually 0.01 is OK, the smaller the more accurate)
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
    :param _e_ph_integ: integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from _mesh
    :param _rand_meth: method for generation of pseudo-random numbers for e-beam phase-space integration:
        1- standard pseudo-random number generator
        2- Halton sequences
        3- LPtau sequences (to be implemented)
    """

    nProc = 1
    rank = 1
    MPI = None
    comMPI = None
    try:
        #DEBUG
        #resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], -1) #MPI module load
        resImpMPI4Py = __import__('mpi4py', globals(), locals(), ['MPI'], 0) #MPI module load
        #print('__import__ passed')

        MPI = resImpMPI4Py.MPI
        comMPI = MPI.COMM_WORLD
        rank = comMPI.Get_rank()
        nProc = comMPI.Get_size()

    except:
        print('Calculation will be sequential (non-parallel), because "mpi4py" module can not be loaded')

    #print(MPI)
    #print(rank, nProc)

    #if(nProc <= 1): #OC050214
    #    _n_part_avg_proc = _n_part_tot

    wfr = SRWLWfr() #Wavefronts to be used in each process
    wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
    wfr.mesh.set_from_other(_mesh)    
    wfr.partBeam = deepcopy(_e_beam)
    arPrecParSR = [_sr_meth, _sr_rel_prec, 0, 0, 50000, 0, _sr_samp_fact]

    #meshRes = SRWLRadMesh()
    meshRes = SRWLRadMesh(_mesh.eStart, _mesh.eFin, _mesh.ne, _mesh.xStart, _mesh.xFin, _mesh.nx, _mesh.yStart, _mesh.yFin, _mesh.ny, _mesh.zStart) #to ensure correct final mesh if _opt_bl==None

    ePhIntegMult = 1
    if(_e_ph_integ == 1): #Integrate over Photon Energy
        eAvg = 0.5*(_mesh.eStart + _mesh.eFin)
        ePhIntegMult = 1000*(_mesh.eFin - _mesh.eStart)/eAvg #To obtain photon energy integrated Intensity in [ph/s/mm^2] assuming monochromatic Spectral Intensity in [ph/s/.1%bw/mm^2]
        wfr.mesh.eStart = eAvg
        wfr.mesh.eFin = eAvg
        wfr.mesh.ne = 1
        meshRes.eStart = eAvg
        meshRes.eFin = eAvg
        meshRes.ne = 1

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

    _sr_rel_prec = int(_sr_rel_prec)
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
        nPartPerProc = _n_part_avg_proc*nSentPerProc #Number of electrons treated by each worker process

    useGsnBmSrc = False
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

    resStokes = None
    workStokes = None
    iAvgProc = 0
    iSave = 0

    doMutual = 0
    if((_char >= 2) and (_char <= 4)): doMutual = 1
    
    if(((rank == 0) or (nProc == 1)) and (_opt_bl != None)): #calculate once the central wavefront in the master process (this has to be done only if propagation is required)

        if(useGsnBmSrc):
            srwl.CalcElecFieldGaussian(wfr, _mag, arPrecParSR)
        else:
            srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR)
            
        #print('DEBUG MESSAGE: Central Wavefront calculated')
        srwl.PropagElecField(wfr, _opt_bl)
        #print('DEBUG MESSAGE: Central Wavefront propagated')
        if(_pres_ang > 0):
            srwl.SetRepresElecField(wfr, 'a')
        
        meshRes.set_from_other(wfr.mesh)

        if(doMutual > 0):
            if(_char == 2):
                meshRes.ny = 1
                meshRes.yStart = _y0
                meshRes.yFin = _y0
            elif(_char == 3):
                meshRes.nx = 1
                meshRes.xStart = _x0
                meshRes.xFin = _x0

        if(nProc > 1): #send resulting mesh to all workers
            #comMPI.send(wfr.mesh, dest=)
            arMesh = array('f', [meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny])
            #comMPI.Bcast([arMesh, MPI.FLOAT], root=MPI.ROOT)
            #comMPI.Bcast([arMesh, MPI.FLOAT])

            for iRank in range(nProc - 1):
                dst = iRank + 1
                #print("msg %d: sending data from %d to %d" % (iRank, rank, dst)) #an he
                comMPI.Send([arMesh, MPI.FLOAT], dest=dst)
            #print('DEBUG MESSAGE: Mesh of Propagated central wavefront broadcasted')

        #DEBUG
        #print('meshRes: ne=', meshRes.ne, 'eStart=', meshRes.eStart, 'eFin=', meshRes.eFin)
        #END DEBUG

        resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
        #wfr.calc_stokes(resStokes) #OC190414 (don't take into account first "central" beam)
        workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

        #iAvgProc += 1 #OC190414 (commented-out)
        #iSave += 1

    #slaves = [] #an he
    #print('DEBUG MESSAGE: rank=', rank)
    if((rank > 0) or (nProc == 1)):

        if((nProc > 1) and (_opt_bl != None)): #receive mesh for the resulting wavefront from the master
            arMesh = array('f', [0]*9)
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
            #sys.exit(0)

        nRadPt = meshRes.ne*meshRes.nx*meshRes.ny
        if(doMutual > 0): nRadPt *= nRadPt
        
        nStPt = nRadPt*4
        randAr = array('d', [0]*6) #for random Gaussian numbers

        #random.seed(rank)
        random.seed(rank*123)
        newSeed = random.randint(0, 1000000)
        random.seed(newSeed)
        
        iAuxSendCount = 0 #for debug
        
        for i in range(nPartPerProc): #loop over macro-electrons

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
            wfr.partBeam.partStatMom1.gamma = elecGamma0*(1 + elecAbsEnSpr*randAr[4]/elecE0)

            #reset mesh, because it may be modified by CalcElecFieldSR and PropagElecField
            wfr.mesh.set_from_other(_mesh)

            if(_e_ph_integ == 1):
                if(_rand_meth == 1):
                    ePh = random.uniform(_mesh.eStart, _mesh.eFin)
                else:
                    ePh = _mesh.eStart + (_mesh.eFin - _mesh.eStart)*randAr[5]
                    
                wfr.mesh.eStart = ePh
                wfr.mesh.eFin = ePh
                wfr.mesh.ne = 1
            
            wfr.presCA = 0 #presentation/domain: 0- coordinates, 1- angles
            wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time

            if(nProc == 1):
                print('i=', i, 'Electron Coord.: x=', wfr.partBeam.partStatMom1.x, 'x\'=', wfr.partBeam.partStatMom1.xp, 'y=', wfr.partBeam.partStatMom1.y, 'y\'=', wfr.partBeam.partStatMom1.yp, 'E=',  wfr.partBeam.partStatMom1.gamma*0.51099890221e-03)
                if(_e_ph_integ == 1):
                     print('Eph=', wfr.mesh.eStart)

            try:
                if(useGsnBmSrc):
                    _mag.x = wfr.partBeam.partStatMom1.x
                    _mag.xp = wfr.partBeam.partStatMom1.xp
                    _mag.y = wfr.partBeam.partStatMom1.y
                    _mag.yp = wfr.partBeam.partStatMom1.yp
                    srwl.CalcElecFieldGaussian(wfr, _mag, arPrecParSR)
                    #print('Gaussian wavefront calc. done')
                else:
                    srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecParSR) #calculate Electric Field emitted by current electron

                if(_opt_bl != None):
                    srwl.PropagElecField(wfr, _opt_bl) #propagate Electric Field emitted by the electron

                if(_pres_ang > 0):
                    srwl.SetRepresElecField(wfr, 'a')

            except:
                traceback.print_exc()

            meshWork = deepcopy(wfr.mesh)

            if(doMutual > 0):
                if(_char == 2):
                    meshWork.ny = 1
                    meshWork.yStart = _y0
                    meshWork.yFin = _y0
                elif(_char == 3):
                    meshWork.nx = 1
                    meshWork.xStart = _x0
                    meshWork.xFin = _x0

            if(workStokes == None):
                workStokes = SRWLStokes(1, 'f', meshWork.eStart, meshWork.eFin, meshWork.ne, meshWork.xStart, meshWork.xFin, meshWork.nx, meshWork.yStart, meshWork.yFin, meshWork.ny, doMutual)
            else:
                nRadPtCur = meshWork.ne*meshWork.nx*meshWork.ny
                if(doMutual > 0):
                    nRadPtCur *= nRadPtCur

                nPtCur = nRadPtCur*4
                
                if(len(workStokes.arS) < nPtCur):
                    del workStokes.arS
                    workStokes.arS = array('f', [0]*nPtCur)
                    #workStokes.mesh.set_from_other(wfr.mesh)

            wfr.calc_stokes(workStokes) #calculate Stokes parameters from Electric Field
            #DEBUG
            #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, 1)
            #END DEBUG

            if(resStokes == None):
                resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)
                #DEBUG
                #print('resStokes #2: ne=', resStokes.mesh.ne, 'eStart=', resStokes.mesh.eStart, 'eFin=', resStokes.mesh.eFin)
                #END DEBUG
 
            if(_opt_bl == None):
                #resStokes.avg_update_same_mesh(workStokes, iAvgProc, 1)
                resStokes.avg_update_same_mesh(workStokes, iAvgProc, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
                
                #DEBUG
                #srwl_uti_save_intens_ascii(workStokes.arS, workStokes.mesh, _file_path, 1)
                #END DEBUG
                
            else:
                #print('DEBUG MESSAGE: Started interpolation of current wavefront on resulting mesh')
                #if(doMutual <= 0): resStokes.avg_update_interp(workStokes, iAvgProc, 1, 1)
                #else: resStokes.avg_update_interp_mutual(workStokes, iAvgProc, 1)

                if(doMutual <= 0): resStokes.avg_update_interp(workStokes, iAvgProc, 1, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
                else: resStokes.avg_update_interp_mutual(workStokes, iAvgProc, 1, ePhIntegMult)
                
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

                    comMPI.Send([resStokes.arS, MPI.FLOAT], dest=0)
                    iAuxSendCount += 1 #for debug

                    #DEBUG
                    #srwl_uti_save_text("Sent # " + str(iAuxSendCount), _file_path + "." + str(rank) + "es.dbg")
                    #END DEBUG

                    for ir in range(nStPt):
                        resStokes.arS[ir] = 0
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
                if((_file_path != None) and (iSave == _n_save_per)):
                    #Saving results from time to time in the process of calculation:
                    srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
                    #sys.exit(0)
                    iSave = 0

    elif((rank == 0) and (nProc > 1)):

        #nRecv = int(nPartPerProc*nProc/_n_part_avg_proc + 1e-09)
        nRecv = nSentPerProc*(nProc - 1) #Total number of sending acts to be made by all worker processes, and to be received by master
        
        #DEBUG
        #srwl_uti_save_text("nRecv: " + str(nRecv) + " nPartPerProc: " + str(nPartPerProc) + " nProc: " + str(nProc) + " _n_part_avg_proc: " + str(_n_part_avg_proc), _file_path + ".00.dbg")
        #END DEBUG

        if(resStokes == None):
            resStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

        if(workStokes == None):
            workStokes = SRWLStokes(1, 'f', meshRes.eStart, meshRes.eFin, meshRes.ne, meshRes.xStart, meshRes.xFin, meshRes.nx, meshRes.yStart, meshRes.yFin, meshRes.ny, doMutual)

        for i in range(nRecv): #loop over messages from workers

            #DEBUG
            #srwl_uti_save_text("Preparing to receiving # " + str(i), _file_path + ".br.dbg")
            #END DEBUG
           
            comMPI.Recv([workStokes.arS, MPI.FLOAT], source=MPI.ANY_SOURCE) #receive #an he (commented-out)

            #DEBUG
            #srwl_uti_save_text("Received intensity # " + str(i), _file_path + ".er.dbg")
            #END DEBUG

            #resStokes.avg_update_same_mesh(workStokes, i + 1)
            #resStokes.avg_update_same_mesh(workStokes, i + 1, 1, ePhIntegMult) #to treat all Stokes components / Polarization in the future
            multFinAvg = 1 if(_n_part_avg_proc > 1) else ePhIntegMult #OC120714 fixed: the normalization may have been already applied at the previous avaraging on each worker node!
            resStokes.avg_update_same_mesh(workStokes, i + 1, 1, multFinAvg) #in the future treat all Stokes components / Polarization, not just s0!

            #DEBUG
            #srwl_uti_save_text("Updated Stokes after receiving intensity # " + str(i), _file_path + "." + str(i) + "er.dbg")
            #END DEBUG

            iSave += 1
            if(iSave == _n_save_per):
                #Saving results from time to time in the process of calculation
                srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)  
                iSave = 0

    #DEBUG
    #srwl_uti_save_text("Exiting srwl_wfr_emit_prop_multi_e", _file_path + "." + str(rank) + "e.dbg")
    #END DEBUG

    if((rank == 0) or (nProc == 1)):
        #Saving final results:
        if(_file_path != None):
            srwl_uti_save_intens_ascii(resStokes.arS, meshRes, _file_path, 1, _mutual = doMutual)
            
        return resStokes
    else:
        return None

#****************************************************************************
#****************************************************************************
#Import of modules requiring classes defined in this smodule
#****************************************************************************
#****************************************************************************
from srwl_uti_src import *

