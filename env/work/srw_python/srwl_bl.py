# -*- coding: utf-8 -*-
#############################################################################
# SRWLib SR Beamline Base Class
# Contains a set of member objects and functions for simulating basic operation and characteristics
# of a complete user beamline in a synchrotron radiation source.
# Under development!!!
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import uti_math
import time

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
    def set_e_beam(self, _e_beam_name='', _e_beam=None, _i=-1, _sig_e=-1, _emit_x=-1, _emit_y=-1, _drift=0, _x=0, _y=0, _xp=0, _yp=0, _dE=0):
        """Setup Electron Beam
        :param _e_beam_name: e-beam unique name, e.g. 'NSLS-II Low Beta Day 1' (see srwl_uti_src.py)
        :param _e_beam: e-beam structure (SRWLPartBeam instance)
        :param _i: e-beam current [A]
        :param _sig_e: e-beam relative energy spread
        :param _emit_x: e-beam horizontal emittance
        :param _emit_y: e-beam vertical emittance
        :param _drift: e-beam drift length in [m] from center of straight section
        :param _x: initial average horizontal position [m]
        :param _y: initial average vertical position [m]
        :param _xp: initial average horizontal angle [m]
        :param _yp: initial average vertical angle [m]
        :param _dE0: average energy deviation [GeV]
        """
        #add more parameters when/if necessary

        if(_sig_e < 0.): _sig_e = None
        if(_emit_x < 0.): _emit_x = None
        if(_emit_y < 0.): _emit_y = None

        sIncInpElecBeam = 'Incorrect input for setting up Electron Beam structure'
        if(len(_e_beam_name) > 0):
            self.eBeam = srwl_uti_src_e_beam(_e_beam_name, _sig_e=_sig_e, _emit_x=_emit_x, _emit_y=_emit_y)
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
        self.eBeam.partStatMom1.xp = _xp
        self.eBeam.partStatMom1.yp = _yp
        
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

        #print(self.mag_approx.arVx)
        
        return self.mag_approx

    #------------------------------------------------------------------------
    def set_und_tab(self, _gap, _ph_mode='p1', _phase=0., _zc=0., _interp_ord=1, _meas_or_calc='m'):
        """Setup magnetic field container with magnetic measurements or calculation data interpolated for given gap and phase
        :param _gap: magnetic gap [mm] for which the field should be set up
        :param _ph_mode: type of phase (shift of magnet arrays) motion
        :param _phase: shift of magnet arrays [mm] for which the field should be set up
        :param _zc: center position [m]
        :param _interp_ord: order of interpolation: 1- (bi-)linear, 2- (bi-)quadratic, 3- (bi-)cubic
        :param _meas_or_calc: use magnetic measurements ('m') or calculation ('c') data
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
        
        for i in range(nRows):
            curLine = lines[i]
            curLineParts = curLine.split(strSep)
            curLenLineParts = len(curLineParts)
            if(curLenLineParts >= 4):
                curPhaseMode = curLineParts[1]
                if(curPhaseMode != _ph_mode): continue

                arGaps.append(float(curLineParts[0]))

                curPhase = float(curLineParts[2])
                if((phasePrev != None) and (curPhase != phasePrev)): phaseIsVar = True
                arPhases.append(curPhase)
                phasePrev = curPhase
                
                curFileName = curLineParts[3]
                #print(curFileName)
                
                if(len(curFileName) > 0):
                    curFilePath = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, curFileName)
                    curFldCnt = srwl_uti_read_mag_fld_3d(curFilePath, '#')
                    arMagFld3D.append(curFldCnt.arMagFld[0])
                    arXc.append(curFldCnt.arXc[0])
                    arYc.append(curFldCnt.arYc[0])
                    arZc.append(curFldCnt.arZc[0] + _zc)
                if(curLenLineParts >= 6):
                    arCoefBx.append(float(curLineParts[4]))
                    arCoefBy.append(float(curLineParts[5]))
        f.close()

        fldCnt = SRWLMagFldC(arMagFld3D, array('d', arXc), array('d', arYc), array('d', arZc))
        nElem = len(arMagFld3D)
        if((nElem != len(arGaps)) or (nElem != len(arPhases))):
            raise Exception('Inconsistent magnetic field data summary file')
        
        fldCnt.arPar1 = array('d', arGaps)
        fldCnt.arPar2 = array('d', arPhases)
        if(len(arCoefBx) == nElem): fldCnt.arPar3 = array('d', arCoefBx)
        if(len(arCoefBy) == nElem): fldCnt.arPar4 = array('d', arCoefBy)

        fldCntRes = SRWLMagFldC(arMagFld3D[0], arXc[0], arYc[0], arZc[0])
        precPar = [1, _gap, _phase, _interp_ord]
        self.mag = srwl.CalcMagnField(fldCntRes, fldCnt, precPar)
        return self.mag

    #------------------------------------------------------------------------
    def set_und_per_from_tab(self, _rel_ac_thr=0.05, _max_nh=5, _max_per=0.1):
        """Setup periodic Magnetic Field from Tabulated one
        :param _rel_ac_thr: relative accuracy threshold
        :param _max_nh: max. number of harmonics to create
        :param _max_per: max. period length to consider
        """

        sErMes = 'Magnetic Field is not defined'
        if(self.mag == None): Exception(sErMes)
        if(isinstance(self.mag, SRWLMagFldC) == False): raise Exception(sErMes)

        arHarm = []
        for i in range(_max_nh): arHarm.append(SRWLMagFldH())

        self.mag_approx = SRWLMagFldC(SRWLMagFldU(arHarm))
        srwl.UtiUndFromMagFldTab(self.mag_approx, self.mag, [_rel_ac_thr, _max_nh, _max_per])
        return self.mag_approx

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

        if(self.gsnBeam != None): del self.gsnBeam

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
    def set_optics(self, _op):
        """Setup optical element container
        :param _op: optical element container (SRWLOptC instance)
        """

        if((_op == None) or (isinstance(_op, SRWLOptC) == False)):
            raise Exception('Incorrect optics container (SRWLOptC) structure')
        if(self.optics != None): del self.optics
        self.optics = _op

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

        if(self.eBeam == None): Exception('Electron Beam structure is not defined')

        if(_mag_type == 1):
            if(self.mag_approx == None): Exception('Approximate Magnetic Field is not defined')
        elif(_mag_type == 2):
            if(self.mag == None): Exception('Magnetic Field is not defined')
        else: Exception('Incorrect Magnetic Field type identificator')

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

        print('Electron trajectory calculatiton ... ', end='')
        srwl.CalcPartTraj(partTraj, magToUse, arPrecPar)
        print('completed')

        if(len(_fname) > 0):
            print('Saving trajectory data to a file ... ', end='')
            partTraj.save_ascii(_fname)
            print('completed')

        return partTraj

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
        if(_mag_type == 2):
            magToUse = self.mag
            #print('Using tabulated magnetic field...')

        wfr = SRWLWfr()
        wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh = deepcopy(_mesh)
        wfr.partBeam = self.eBeam

        zStartInteg = 0 #longitudinal position to start integration (effective if < zEndInteg)
        zEndInteg = 0 #longitudinal position to finish integration (effective if > zStartInteg)
        npTraj = 50000 #Number of points for trajectory calculation 
        useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
        arPrecPar = [_meth, _rel_prec, zStartInteg, zEndInteg, npTraj, useTermin, _samp_fact]

        #print('magToUse=', magToUse.arMagFld[0])

        print('Single-electron SR calculatiton ... ', end='')
        t0 = time.time();
        srwl.CalcElecFieldSR(wfr, 0, magToUse, arPrecPar) #calculate SR
        print('completed (lasted', round(time.time() - t0, 6), 's)')

        arI = None
        if(_int_type >= 0):
            print('Extracting intensity and saving it to a file ... ', end='')
            t0 = time.time();
            sNumTypeInt = 'f'
            if(_int_type == 4): sNumTypeInt = 'd'

            arI = array(sNumTypeInt, [0]*wfr.mesh.ne*wfr.mesh.nx*wfr.mesh.ny)
            srwl.CalcIntFromElecField(arI, wfr, _pol, _int_type, depType, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
            if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, wfr.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
            print('completed (lasted', round(time.time() - t0, 6), 's)')
        return wfr, arI

    #------------------------------------------------------------------------
    def calc_rad_gsn(self, _mesh, _samp_fact=-1, _pol=6, _int_type=0, _presFT='f', _unitE=2, _fname=''):
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

        wfr = SRWLWfr()
        wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
        wfr.mesh = deepcopy(_mesh)

        wfr.presFT = 0 #presentation/domain: 0- frequency (photon energy), 1- time
        if(_presFT == "t"): wfr.presFT = 1

        wfr.unitElFld = _unitE;

        wfr.partBeam.partStatMom1.x = self.gsnBeam.x #Some information about the source in the Wavefront structure
        wfr.partBeam.partStatMom1.y = self.gsnBeam.y
        wfr.partBeam.partStatMom1.z = self.gsnBeam.z
        wfr.partBeam.partStatMom1.xp = self.gsnBeam.xp
        wfr.partBeam.partStatMom1.yp = self.gsnBeam.yp

        print('Gaussian beam electric field calculatiton ... ', end='')
        t0 = time.time();
        srwl.CalcElecFieldGaussian(wfr, self.gsnBeam, [_samp_fact])
        print('completed (lasted', round(time.time() - t0, 6), 's)')

        arI = None
        if(_int_type >= 0):
            print('Extracting intensity and saving it to a file ... ', end='')
            t0 = time.time();
            sNumTypeInt = 'f'
            if(_int_type == 4): sNumTypeInt = 'd'

            #print('depType=', depType)

            arI = array(sNumTypeInt, [0]*wfr.mesh.ne*wfr.mesh.nx*wfr.mesh.ny)
            srwl.CalcIntFromElecField(arI, wfr, _pol, _int_type, depType, wfr.mesh.eStart, wfr.mesh.xStart, wfr.mesh.yStart)
            if(len(_fname) > 0): srwl_uti_save_intens_ascii(arI, wfr.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', 'ph/s/.1%bw/mm^2'])
            print('completed (lasted', round(time.time() - t0, 6), 's)')
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
    def calc_arb_spec_me(self, _mesh, _meth=2, _rel_prec=0.01, _n_part_tot=100000, _n_part_avg_proc=10, _n_save_per=10, _type=2, _mag=2, _pol=0, _rand_meth=1, _fname=None):
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
        :return: 1D array with (C-aligned) resulting intensity data
        """

        if((_mesh == None) or (isinstance(_mesh, SRWLRadMesh) == False)):
            raise Exception('Incorrect SRWLRadMesh structure')

        if(self.eBeam == None): Exception('Electron Beam structure is not defined')

        mag2use = self.mag
        if(_mag == 1):
            if(self.mag_approx == None): Exception('Approximate Magnetic Field is not defined')
            mag2use = self.mag_approx
        else:
            if(self.mag == None): Exception('Accurate Magnetic Field is not defined')

        charMultiE = 0 #Calculate intensity (flux per unit surface by default)
        if(_type == 1): charMultiE = 10 #Calculate flux

        stk = srwl_wfr_emit_prop_multi_e(
            _e_beam = self.eBeam, _mag = mag2use, _mesh = _mesh,
            _sr_meth = _meth, _sr_rel_prec = _rel_prec,
            _n_part_tot = _n_part_tot, _n_part_avg_proc = _n_part_avg_proc, _n_save_per = _n_save_per, _rand_meth = _rand_meth,
            _file_path = _fname, _char = charMultiE)

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
    def calc_wfr_prop(self, _wfr, _pres_ang=0, _pol=6, _int_type=0, _dep_type=3, _fname=''):
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
        :return: 1D array with (C-aligned) resulting intensity data; it also modified _wfr in place
        """

        if((hasattr(self, 'optics') == False) or (isinstance(self.optics, SRWLOptC) == False)):
            raise Exception('Incorrect optics container (SRWLOptC) structure')
        if(isinstance(_wfr, SRWLWfr) == False):
            raise Exception('Incorrect wavefront (SRWLWfr) structure')

        print('Propagation ... ', end='')
        t0 = time.time();
        srwl.PropagElecField(_wfr, self.optics)
        print('completed (lasted', round(time.time() - t0, 6), 's)')

        if(_pres_ang != 0): srwl.SetRepresElecField(_wfr, 'a')
 
        arI = None
        if(_int_type >= 0): 
            sNumTypeInt = 'f'
            if(_int_type == 4): sNumTypeInt = 'd'

            arI = array(sNumTypeInt, [0]*_wfr.mesh.ne*_wfr.mesh.nx*_wfr.mesh.ny)
            srwl.CalcIntFromElecField(arI, _wfr, _pol, _int_type, _dep_type, _wfr.mesh.eStart, _wfr.mesh.xStart, _wfr.mesh.yStart)
            if(len(_fname) > 0):
                sValUnitName = 'ph/s/.1%bw/mm^2' #consider allowing for other units (for FEL applications)

                print('Saving Propagation Results ... ', end='')
                t0 = time.time();
                srwl_uti_save_intens_ascii(arI, _wfr.mesh, _fname, 0, ['Photon Energy', 'Horizontal Position', 'Vertical Position', ''], _arUnits=['eV', 'm', 'm', sValUnitName])
                print('completed (lasted', round(time.time() - t0), 's)')

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

        #---main folder
        if hasattr(_v, 'fdir'): self.dir_main = _v.fdir

        #---setup electron beam
        if(hasattr(_v, 'ebm_nm')): #To improve
            self.set_e_beam(
                _e_beam_name = (_v.ebm_nm + _v.ebm_nms),
                _i = _v.ebm_i,
                _sig_e = _v.ebm_ens,
                _emit_x = _v.ebm_emx,
                _emit_y = _v.ebm_emy,
                _drift = _v.ebm_dr,
                _x = _v.ebm_x,
                _y = _v.ebm_y,
                _xp = _v.ebm_xp,
                _yp = _v.ebm_yp,
                _dE = _v.ebm_de)

        #print('e-beam was set up')

        #---setup magnetic field: undulator, sinusoidal approximation
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
            self.mag = None
            
            if _v.und_b2e:
                e1 = self.mag_approx.arMagFld[len(self.mag_approx.arMagFld) - 1].get_E1(_en_elec=self.eBeam.partStatMom1.get_E(), _unit='eV')
                print('Fundamental Photon Energy:', srwl_uti_num_round(e1), 'eV') #check how it will work under IPython

            if _v.und_e2b:
                b = self.mag_approx.arMagFld[len(self.mag_approx.arMagFld) - 1].E1_2_B(_e1=_v.w_e, _en_elec=self.eBeam.partStatMom1.get_E())
                print('Magnetic Field Amplitude:', srwl_uti_num_round(b), 'T') #check how it will work under IPython

        #---setup magnetic field: undulator, tabulated magnetic field (measured)
        magnMeasDirExists = False
        if hasattr(_v, 'und_mdir'):
            self.dir_magn_meas = _v.und_mdir
            if(len(self.dir_magn_meas) > 0): magnMeasDirExists = True
        magnMeasSumFileExists = False
        if hasattr(_v, 'und_mfs'):
            self.fn_magn_meas_sum = _v.und_mfs
            if(magnMeasDirExists):
                testPath = os.path.join(os.getcwd(), self.dir_main, self.dir_magn_meas, self.fn_magn_meas_sum)
                magnMeasSumFileExists = os.path.exists(testPath) 

        if magnMeasSumFileExists and hasattr(_v, 'und_g'):
            if(_v.und_g > 0.):
                phase = 0.
                if hasattr(_v, 'und_ph'): phase = _v.und_ph
                phase_mode = 'p1'
                if hasattr(_v, 'und_phm'): phase_mode = _v.und_phm

                self.set_und_tab(#setup undulator source from measured magnetic field data
                    _gap = _v.und_g,
                    _ph_mode = phase_mode,
                    _phase = phase,
                    _zc = _v.und_zc,
                    _interp_ord = 3, #1,
                    _meas_or_calc='m')

                #if((_v.ss_mag == 1) or (_v.ss_mag == 1) or (_v.w_mag == 1) or (_v.tr_mag == 1)):
                if((_v.ss and (_v.ss_mag == 1)) or (_v.sm and (_v.sm_mag == 1)) or ((_v.ws or _v.wm) and (_v.w_mag == 1)) or (_v.tr and (_v.tr_mag == 1))):
                    #print('test')
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

        #---setup Gaussian beam
        if hasattr(_v, 'gbm_pen'):
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

        #---calculate electron trajectory
        if(_v.tr): 
            trj = self.calc_el_trj(
                _ctst = _v.tr_cti, _ctfi = _v.tr_ctf, _np = _v.tr_np,
                _mag_type = _v.tr_mag,
                _fname = os.path.join(_v.fdir, _v.tr_fn) if(len(_v.tr_fn) > 0) else '')

        #---calculate single-e spectrum vs photon energy
        if(_v.ss or _v.gs): 
            mesh_ss = SRWLRadMesh(
                _v.ss_ei, _v.ss_ef, _v.ss_ne,
                _v.ss_x, _v.ss_x, 1,
                _v.ss_y, _v.ss_y, 1,
                _v.op_r)
            
            srCanBeCalc = (self.eBeam != None) and ((self.mag_approx != None) or (self.mag != None))
            #print("                    _v.ss=", _v.ss)
            #print("                    _v.gs=", _v.gs)
            
            if((_v.gs != True) and (srCanBeCalc == True)):
                #print("                                  Before calc_sr_se")
                #print("                                  self.mag=", self.mag)
                wfr_ss, int_ss = self.calc_sr_se(
                    _mesh = mesh_ss,
                    _meth = _v.ss_meth,
                    _rel_prec = _v.ss_prec,
                    _pol = _v.ss_pol,
                    _int_type = 0,
                    _mag_type = _v.ss_mag,
                    _fname = os.path.join(_v.fdir, _v.ss_fn) if(len(_v.ss_fn) > 0) else '')
            if((_v.gs == True) or ((self.gsnBeam != None) and (srCanBeCalc == False))):
                wfr_ss, int_ss = self.calc_rad_gsn(
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
                int_sm = self.calc_ur_spec_me(
                    _mesh = mesh_sm,
                    _harm_init = _v.sm_hi,
                    _harm_fin = _v.sm_hf,
                    _prec_long = _v.sm_prl,
                    _prec_azim = _v.sm_pra,
                    _type = _v.sm_type,
                    _pol = _v.sm_pol,
                    _fname = os.path.join(_v.fdir, _v.sm_fn) if(len(_v.sm_fn) > 0) else '')
            else: #for "accurate" magnetic field
                int_sm = self.calc_arb_spec_me(
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
                    _fname = os.path.join(_v.fdir, _v.sm_fn) if(len(_v.sm_fn) > 0) else '')

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
                #_z_start = _v.pw_zst,
                #_z_fin = _v.pw_zfi,
                _z_start = _v.pw_zi,
                _z_fin = _v.pw_zf,
                _mag_type = _v.pw_mag,
                _fname= os.path.join(_v.fdir, _v.pw_fn) if(len(_v.pw_fn) > 0) else '')
            
        #---calculate single-e and multi-e intensity distributions (before and after wavefront propagation through a beamline)
        if(_v.si or _v.ws or _v.gi or _v.wg or _v.wm):
            if(_v.ws or _v.wg or _v.wm): self.set_optics(_op)
            
            ef = _v.w_e
            if(_v.w_ef > 0): ef = _v.w_ef
            mesh_w = SRWLRadMesh(
                _v.w_e, ef, _v.w_ne,
                _v.w_x - 0.5*_v.w_rx, _v.w_x + 0.5*_v.w_rx, _v.w_nx,
                _v.w_y - 0.5*_v.w_ry, _v.w_y + 0.5*_v.w_ry, _v.w_ny,
                _v.op_r)
        #---calculate single-e electric field and intensity (before wavefront propagation through a beamline)
            if(_v.si or _v.ws or _v.gi or _v.wg):

                srCanBeCalc = (self.eBeam != None) and ((self.mag_approx != None) or (self.mag != None))
                gsnBeamCanBeCalc = self.gsnBeam != None

                #if((_v.gi == False) and (_v.wg == False) and (srCanBeCalc == True)):
                if((_v.gi != True) and (_v.wg != True) and (srCanBeCalc == True)):
                    wfr, int_w0 = self.calc_sr_se(
                        _mesh = deepcopy(mesh_w),
                        _samp_fact = _v.w_smpf,
                        _meth = _v.w_meth,
                        _rel_prec = _v.w_prec,
                        _pol = _v.si_pol,
                        _int_type = _v.si_type,
                        _mag_type = _v.w_mag,
                        _fname = os.path.join(_v.fdir, _v.si_fn) if(len(_v.si_fn) > 0) else '')

                if((_v.gs == True) or ((gsnBeamCanBeCalc == True) and (srCanBeCalc == False))):
                    wfr, int_w0 = self.calc_rad_gsn(
                        _mesh = deepcopy(mesh_w),
                        _samp_fact = _v.w_smpf,
                        _pol = _v.si_pol,
                        _int_type = _v.si_type,
                        _presFT = _v.w_ft,
                        _unitE = _v.w_u,
                        _fname = os.path.join(_v.fdir, _v.si_fn) if(len(_v.si_fn) > 0) else '')
                    
                mesh_si = deepcopy(wfr.mesh)
                
        #---calculate single-e electric field and intensity (after wavefront propagation through a beamline)
                if(_v.ws or _v.wg):
                    int_ws = self.calc_wfr_prop(
                        _wfr = wfr,
                        _pres_ang = _v.ws_ap,
                        _pol = _v.si_pol,
                        _int_type = _v.si_type,
                        _dep_type=3, #consider adding other cases (e.g. for TD FEL calculations)
                        _fname = os.path.join(_v.fdir, _v.ws_fni) if(len(_v.ws_fni) > 0) else '')
                    mesh_ws = wfr.mesh
                    #if(len(_v.ws_fn) > 0): to implement saving single-e (/ fully coherent) wavefront data (wfr) to a file

        #---calculate multi-electron (/ partially coherent) wavefront propagation
            if(_v.wm):
                wmResFileName = os.path.join(_v.fdir, _v.wm_fni) if(len(_v.wm_fni) > 0) else None
                #print(wmResFileName)
                
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

        if (_v.tr == True) and (len(_v.tr_pl) > 0):
            #uti_plot1d(int_ss, [mesh_ss.eStart, mesh_ss.eFin, mesh_ss.ne], ['Photon Energy', 'Intensity', 'Intensity'], ['eV', 'ph/s/.1%bw/mm^2'])

            #To implement!
            #ddddddddddddddddddddddddddddddddddddddddddddddddd

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
            
            uti_plot1d(int_ss, [mesh_ss.eStart, mesh_ss.eFin, mesh_ss.ne], [sArgLabel, sValLabel, sValLabel], [sArgUnit, sValUnit])
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
                    int_w0,
                    [mesh_si.xStart, mesh_si.xFin, mesh_si.nx],
                    [mesh_si.yStart, mesh_si.yFin, mesh_si.ny],
                    0, #0.5*(mesh_si.xStart + mesh_si.xFin),
                    0, #0.5*(mesh_si.yStart + mesh_si.yFin),
                    ['Horizontal Position', 'Vertical Position', sValLabel],
                    ['m', 'm', sValUnit], #add other units for FEL
                    True)
                plotOK = True

        if _v.ws and (len(_v.ws_pl) > 0):
            if (_v.ws_pl == 'xy') or (_v.ws_pl == 'yx') or (_v.ws_pl == 'XY') or (_v.ws_pl == 'YX'):

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
                    int_w0,
                    [mesh_si.xStart, mesh_si.xFin, mesh_si.nx],
                    [mesh_si.yStart, mesh_si.yFin, mesh_si.ny],
                    0, #0.5*(mesh_si.xStart + mesh_si.xFin),
                    0, #0.5*(mesh_si.yStart + mesh_si.yFin),
                    ['Horizontal Position', 'Vertical Position', sValLabel + ' Before Propagation'],
                    ['m', 'm', sValUnit],
                    True)
                uti_plot2d1d(
                    int_ws,
                    [mesh_ws.xStart, mesh_ws.xFin, mesh_ws.nx],
                    [mesh_ws.yStart, mesh_ws.yFin, mesh_ws.ny],
                    0, #0.5*(mesh_ws.xStart + mesh_ws.xFin),
                    0, #0.5*(mesh_ws.yStart + mesh_ws.yFin),
                    ['Horizontal Position', 'Vertical Position', sValLabel + ' After Propagation'],
                    ['m', 'm', sValUnit],
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
def srwl_uti_std_options():
    """Defines set of standard options (applicable to any beamline) for general types of calculation
    :returns: list providing compact description of all options; every element of this list is supposed to contain:
        [0]: string containing option (/ variable) name
        [1]: string containing type of the option / variable ('f' - float, 'i' - integer, 's' - string)
        [2]: default value
        [3]: string containing help / explanation of the option / variable
        [4]: optional string describing formal action to be taken if option is fired
    """
    varParamStd = [
#---Undulator
        ['und_mdir', 's', 'magn_meas', 'name of magnetic measurements sub-folder'],
        ['und_g', 'f', 0., 'undulator gap [mm] (assumes availability of magnetic measurement or simulation data)'],
        ['und_b2e', '', '', 'estimate undulator fundamental photon energy (in [eV]) for the amplitude of sinusoidal magnetic field defined by und_b or und_bx, und_by', 'store_true'],
        ['und_e2b', '', '', 'estimate undulator field amplitude (in [T]) for the photon energy defined by w_e', 'store_true'],

#---Calculation Types
    #Electron Trajectory
        ['tr', '', '', 'calculate electron trajectory', 'store_true'],
        ['tr_cti', 'f', 0., 'initial time moment (c*t) for electron trajectory calculation [m]'],
        ['tr_ctf', 'f', 0., 'final time moment (c*t) for electron trajectory calculation [m]'],
        ['tr_np', 'f', 50000, 'number of points for trajectory calculation'],
        ['tr_mag', 'i', 2, 'magnetic field to be used for trajectory calculation: 1- approximate, 2- accurate'],
        ['tr_fn', 's', 'res_trj.dat', 'file name for saving calculated trajectory data'],
        ['tr_pl', 's', 'xxpyypz', 'plot the resulting trajectiry in graph(s): ""- dont plot, otherwise the string should list the trajectory components to plot'],

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
        ['sm_pol', 'i', 6, 'polarization component to extract after calculation of multi-e flux or intensity: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
        ['sm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
        ['sm_fn', 's', 'res_spec_me.dat', 'file name for saving calculated milti-e spectrum vs photon energy'],
        ['sm_pl', 's', 'e', 'plot the resulting spectrum-e spectrum in a graph: ""- dont plot, "e"- show plot vs photon energy'],

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
        ['pw_zi', 'f', 0., 'initial longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
        ['pw_zf', 'f', 0., 'final longitudinal position along electron trajectory of power density distribution (effective if pow_sst < pow_sfi)'],
        ['pw_mag', 'i', 1, 'magnetic field to be used for power density calculation: 1- approximate, 2- accurate'],
        ['pw_fn', 's', 'res_pow.dat', 'file name for saving calculated power density distribution'],
        ['pw_pl', 's', 'xy', 'plot the resulting power density distribution in a graph: ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

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
        ['w_meth', 'i', 1, 'method to use for calculation of intensity distribution vs horizontal and vertical position'],
        ['w_prec', 'f', 0.01, 'relative precision for calculation of intensity distribution vs horizontal and vertical position'],
        ['w_mag', 'i', 1, 'magnetic field to be used for calculation of intensity distribution vs horizontal and vertical position: 1- approximate, 2- accurate'],
        ['w_ft', 's', 'f', 'presentation/domain: "f"- frequency (photon energy), "t"- time'],
        ['w_u', 'i', '1', 'electric field units: 0- arbitrary, 1- sqrt(Phot/s/0.1%bw/mm^2), 2- sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)'],
        
        ['si_pol', 'i', 6, 'polarization component to extract after calculation of intensity distribution: 0- Linear Horizontal, 1- Linear Vertical, 2- Linear 45 degrees, 3- Linear 135 degrees, 4- Circular Right, 5- Circular Left, 6- Total'],
        ['si_type', 'i', 0, 'type of a characteristic to be extracted after calculation of intensity distribution: 0- Single-Electron Intensity, 1- Multi-Electron Intensity, 2- Single-Electron Flux, 3- Multi-Electron Flux, 4- Single-Electron Radiation Phase, 5- Re(E): Real part of Single-Electron Electric Field, 6- Im(E): Imaginary part of Single-Electron Electric Field, 7- Single-Electron Intensity, integrated over Time or Photon Energy'],
        ['si_fn', 's', 'res_int_se.dat', 'file name for saving calculated single-e intensity distribution (without wavefront propagation through a beamline) vs horizontal and vertical position'],
        ['ws_fni', 's', 'res_int_pr_se.dat', 'file name for saving propagated single-e intensity distribution vs horizontal and vertical position'],
        ['ws_pl', 's', 'xy', 'plot the propagated radiaiton intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],
        ['ws_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters (/ Intensity distribution): coordinate (0) or angular (1)'],
        ['si_pl', 's', 'xy', 'plot the input intensity distributions in graph(s): ""- dont plot, "x"- vs horizontal position, "y"- vs vertical position, "xy"- vs horizontal and vertical position'],

        ['wm_nm', 'i', 100000, 'number of macro-electrons (coherent wavefronts) for calculation of multi-electron wavefront propagation'],
        ['wm_na', 'i', 5, 'number of macro-electrons (coherent wavefronts) to average on each node at parallel (MPI-based) calculation of multi-electron wavefront propagation'],
        ['wm_ns', 'i', 5, 'saving periodicity (in terms of macro-electrons / coherent wavefronts) for intermediate intensity at multi-electron wavefront propagation calculation'],
        ['wm_ch', 'i', 0, 'type of a characteristic to be extracted after calculation of multi-electron wavefront propagation: #0- intensity (s0); 1- four Stokes components; 2- mutual intensity cut vs x; 3- mutual intensity cut vs y'],
        ['wm_ap', 'i', 0, 'switch specifying representation of the resulting Stokes parameters: coordinate (0) or angular (1)'],
        ['wm_x0', 'f', 0, 'horizontal center position for mutual intensity cut calculation'],
        ['wm_y0', 'f', 0, 'vertical center position for mutual intensity cut calculation'],
        ['wm_ei', 'i', 0, 'integration over photon energy is required (1) or not (0); if the integration is required, the limits are taken from w_e, w_ef'],
        ['wm_rm', 'i', 1, 'method for generation of pseudo-random numbers for e-beam phase-space integration: 1- standard pseudo-random number generator, 2- Halton sequences, 3- LPtau sequences (to be implemented)'],
        ['wm_fni', 's', 'res_int_pr_me.dat', 'file name for saving propagated multi-e intensity distribution vs horizontal and vertical position'],

        #['ws_fn', 's', '', 'file name for saving single-e (/ fully coherent) wavefront data'],
        #['wm_fn', 's', '', 'file name for saving multi-e (/ partially coherent) wavefront data'],
        #to add options
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
def _optparse(_descr, use_sys_argv=True, args=None):  # MR26022016, MR04032016
    """Set and parse command-prompt options from a compact description provided in _descr using optparse (deprecated since Python 2.7).
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

    # MR07032016:
    if use_sys_argv:
        v, args = p.parse_args()  # MR07032016
    else:
        try:
            v, args = p.parse_args(args if args else [])  # MR07032016
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


def _argparse(_descr, use_sys_argv=True, args=None):  # MR26022016, MR04032016
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

    p = argparse.ArgumentParser(None if use_sys_argv else __name__)  # MR07032016
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

        curOpt[3] = curOpt[3].replace('%', '%%')  # screen special '%' symbol
        if (len(sTypeShort) <= 0):
            p.add_argument('--' + curOpt[0], default=defVal, help=curOpt[3], action=sAct)
        else:
            p.add_argument('--' + curOpt[0], type=sType, default=defVal, help=curOpt[3], action=sAct)

    # MR07032016:
    if use_sys_argv:
        v = p.parse_args()  # MR07032016
    else:
        try:
            v = p.parse_args(args if args else [])  # MR07032016
        except SystemExit as e:
            raise ValueError('Exit code: {}'.format(e))

    # "post-parsing" list-type options
    for i in range(len(listOptNamesPostParse)):
        curOptName = listOptNamesPostParse[i]
        valCurOpt = getattr(v, curOptName)

        if not isinstance(valCurOpt, list) and not isinstance(valCurOpt, array):  # MR07032016
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

srwl_uti_parse_options = _optparse if os.getenv('SRWL_OPTPARSE') else _argparse  # MR07032016
