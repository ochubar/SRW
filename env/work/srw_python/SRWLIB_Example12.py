# -*- coding: utf-8 -*-
#############################################################################
# SRWLIB Example # 12: Simulating Wavefront Propagation through initial part of a Soft X-Ray Undulator Radiation Beamline containing Variable Line Spacing (VLS) Grating
# Based on input and comtributions of N. Canestrari, E. Vescovo, V. Bisogni (BNL)
# v 0.03
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *
from uti_plot import *
import os
import sys
import time

#*********************************Story
help_str='''SRWLIB Python Example # 12:
Simulating emission and propagation of Soft X-ray (250 eV) undulator radiation (UR) wavefront through a beamline containing VLS grating
!!!!!Under development!!!!!
Depending on input options, either single-electron UR wavefront at a fixed photon energy can be calculated and propagated through
the optical scheme, or the calculation of partially-coherent UR from entire electron beam can be started as a loop over "macro-electrons",
using "srwl_wfr_emit_prop_multi_e" function. 
This function can run either in "normal" sequential mode, or in parallel mode under "mpi4py". 
For this, an MPI2 package and the "mpi4py" Python package have to be installed and configured, and this example has to be started e.g. as:
    mpiexec -n 5 python SRWLIB_Example12.py -m 10000
or with some non-default options as follows:
    mpiexec -n 5 python SRWLIB_Example12.py -m 10000 -e 250 -b -0.003 -w 0.00001
For more information on parallel calculations under "mpi4py" please see documentation to the "mpi4py" and MPI.
Note that the long-lasting partially-coherent UR calculation saves from time to time instant average intensity to an ASCII file, 
so the execution of the long loop over "macro-electrons" can be aborted after some time without the danger that all results will be lost.
be sure that a folder data_esm is present in the current working directory.
'''

#*********************************File Names
strExDataFolderName = 'data_example_12' #example data sub-folder name
strIntSE_OutFilePath = os.path.join(os.getcwd(), strExDataFolderName, 'ex12_res_int_se') #file name for output initial single-electron SR intensity data
strIntPropSE_OutFilePath = os.path.join(os.getcwd(), strExDataFolderName, 'ex12_res_int_prop_se') #file name for output propagated single-electron SR intensity data
strIntPropME_OutFilePath = os.path.join(os.getcwd(), strExDataFolderName, 'ex12_res_int_prop_me') #file name for output propagated multi-electron SR intensity data

#*********************************Undulator	
def setUndulator(_en, _bd=0):
  
  BbE  = { #undulator peak magnetic field by energy
    250. : 0.594024,
    450. : 0.405778,
    1000.: 0.187782,
    1500.: 0.375672, 
    2000.: 0.296978
  }

  numPer = 61.5 #Number of ID Periods (without counting for terminations (3.5m, 57mm)
  undPer = 0.057 #Period Length [m]
  By =  BbE[_en]*(1 + _bd) #Peak Vertical field [T]
  print('By=', By, 'T')
  
  phBy =  0 #Initial Phase of the Vertical field component
  sBy  =  1 #Symmetry of the Vertical field component vs Longitudinal position
  xcID =  0 #Transverse Coordinates of Undulator Center [m]
  ycID =  0
  zcID =  0 #0 #Longitudinal Coordinate of Undulator Center [m]

  und = SRWLMagFldU([SRWLMagFldH(1, 'v', By, phBy, sBy, 1)], undPer, numPer) #Planar Undulator
  mag = SRWLMagFldC([und], array('d', [xcID]), array('d', [ycID]), array('d', [zcID])) #Container of all Field elements
  return mag

#*********************************Electron Beam
def setElecBeam(_dist=0):
  ebeam = SRWLPartBeam()
  ebeam.Iavg = 0.5 #Average Current [A]
  #1st order statistical moments
  ebeam.partStatMom1.x  =  0. #200e-06 #Initial Transverse Coordinates (initial Longitudinal Coordinate will be defined later on) [m]
  ebeam.partStatMom1.y  =  0. #-20.e-06 #0.
  ebeam.partStatMom1.z  =  0. #distance center of straight section to center of the und (which is set to 0)   
  ebeam.partStatMom1.xp =  0. #Initial Relative Transverse Velocities
  ebeam.partStatMom1.yp =  0.
  ebeam.partStatMom1.gamma = 3./0.51099890221e-03 #Relative Energy
  #2nd order statistical moments
  ebeam.arStatMom2[ 0] = (107.086e-06)**2 #<(x-x0)^2> [m^2]
  ebeam.arStatMom2[ 1] = 0 #<(x-x0)*(x'-x'0)> [m]
  ebeam.arStatMom2[ 2] = (5.13604e-06)**2 #<(x'-x'0)^2> 
  ebeam.arStatMom2[ 3] = (5.21536e-06)**2 #<(y-y0)^2> [m^2]
  ebeam.arStatMom2[ 4] = 0 #<(y-y0)*(y'-y'0)> [m]
  ebeam.arStatMom2[ 5] = (1.53393e-06)**2 #<(y'-y'0)^2>
  ebeam.arStatMom2[10] = (0.89e-03)**2 #<(E-E0)^2>/E0^2
  if(_dist != 0): ebeam.drift(_dist)
  return ebeam

#*********************************Mesh for Initial Radiation Wavefront calculation
def setRadMesh(_en):
  mesh = SRWLRadMesh()
  mesh.zStart =  34.366 #Longitudinal Position [m] from Center of Undulator SR has to be calculated
  mesh.ne     =  1
  mesh.eStart =  _en #Initial Photon Energy [eV]
  mesh.eFin   =  mesh.eStart #Final Photon Energy [eV]
  hor_ap      =  6.76e-03 
  ver_ap      =  20.e-03 
  mesh.nx     =  201
  mesh.xStart = -hor_ap*0.5 #-0.00015 #Initial Horizontal Position [m]
  mesh.xFin   =  hor_ap*0.5 # 0.00015 #Final Horizontal Position [m]
  mesh.ny     =  201
  mesh.yStart = -ver_ap*0.5 #-0.00015 #Initial Vertical Position [m]
  mesh.yFin   =  ver_ap*0.5 # 0.00015 #Final Vertical Position [m]
  return mesh

#*********************************Beamline
def setBeamline(_en, _n_opt_el):

  AbE = { #angle by energy (M2, G_alpha, G_beta)
    250. : (85.3317, 88.4755, 82.1879), 
    450. : (86.5319, 88.8813, 84.1826), 
    1000.: (87.6789, 89.2577, 86.1001), 
    1500.: (88.1060, 89.3957, 86.8163), 
    2000.: (88.3603, 89.4775, 87.243) 
  }  

  deg2rad = 3.14159265359/180.
  wavelength = 1.239842e-6/_en
  #More accurate value seems to be "1.23984193"; however, in C part of SRW "1.239842" is used

  #*****Drift M1 -> M2
  M1_M2  = SRWLOptD(20.434) 

  #*****M2 (plane mirror) related
  grazM2 = (90.-AbE[_en][0])*deg2rad
  tM2, sM2 = 0.430, 0.02
  M2A   = SRWLOptA('r', 'a', sM2, tM2*sin(grazM2)) #M2 Aperture

  #*****Drift M2 -> Grating
  M2_G  = SRWLOptD(0.2)

  #*****Grating (VLS plane) related
  vlsG = [1800., 0.08997, 3.004e-6, 9.73e-11] #Polynomial coefficients for VLS Grating Groove Density
  grazG = (90.-AbE[_en][1])*deg2rad
  
  outG = asin(wavelength*vlsG[0]*1.e+03 - cos(grazG)) #Grating Output Angle
  defG = grazG + outG + 1.57079632679 #Grating Deflection Angle

  #Grating Input and Output Arms, Tangential and Sagittal Dimensions:
  pG, qG, tG, sG = 55., 42.63, 0.200, 0.015 

  qGest = cos(outG)**2/(wavelength*vlsG[1]*1.e+06 - sin(grazG)**2/pG)
 
  print('Grating Inc. Angle for Central Energy:', AbE[_en][1], 'deg.')
  print('Grating Exit Angle:', -outG/deg2rad, 'deg.')
  print('Grating Exit Arm (estimated):', qGest, 'm')

  GA = SRWLOptA('r', 'a', sG, tG*sin(grazG)) #Grating Aperture

  #Grating Substrate (plane mirror, deflecting in vertical plane):
  GS = SRWLOptMirPl(_size_tang=tG, _size_sag=sG, _ap_shape='r',
                    _nvx=0, _nvy=cos(grazG), _nvz=-sin(grazG), _tvx=0, _tvy=sin(grazG))
  #The VLS Grating itself:
  G = SRWLOptG(_mirSub=GS, _m=1, _grDen=vlsG[0], _grDen1=vlsG[1], _grDen2=vlsG[2], _grDen3=vlsG[3])

  #*****Drift Grating -> M3
  G_M3 = SRWLOptD(34.63) 

  #*****M3 (elliptical horizontally-focusing mirror) related
  incM3  = 88.75
  grazM3 = (90.-incM3)*deg2rad
  pM3, qM3, tM3, sM3 = 89.63, 8.006, 0.42, 0.02 #M3 Input and Output Arms, Tangential and Sagittal Dimensions

  M3A  = SRWLOptA('r', 'a', tM3*sin(grazM3), sM3) #M3 Aperture

  M3 = SRWLOptMirEl(_p=pM3, _q=qM3, _ang_graz=grazM3, _size_tang=tM3, _size_sag=sM3,
                    _nvx=cos(grazM3), _nvy=0, _nvz=-sin(grazM3), _tvx=-sin(grazM3), _tvy=0)

  #*****Drift M3 -> Secondary Source Aperture (SA)
  M3_SA = SRWLOptD(8.006)

  #*****Secondary Source Aperture (i.e. exit aperture of the monochromator)
  SA = SRWLOptA('r', 'a', 10.e-03, 10.e-06) #dimensions to check/steer

  #*****M4 (ellipsoid of rotation, horizontally- and vertically-focusing mirror) related
  grazM4 = grazM3
  pM4, qM4, tM4, sM4 = 6.01, 0.991, 0.300, 0.050 #M4 Input and Output Arms, Tangential and Sagittal Dimensions

  M4A = SRWLOptA('r', 'a', tM4*sin(grazM4), sM4) #M4 Aperture

  focLenM4 = pM4*qM4/(pM4 + qM4)
  rSagM4 = 2*focLenM4*sin(grazM4)*(1. - 0.0032) #to tune...

  #M4 as Ellipsoid of Rotation (under testing!!!)
  #M4 = SRWLOptMirEl(_p=pM4, _q=qM4, _ang_graz=grazM4, _r_sag=rSagM4, _size_tang=tM4, _size_sag=sM4,
  #                  _nvx=cos(grazM4), _nvy=0, _nvz=-sin(grazM4), _tvx=-sin(grazM4), _tvy=0) 
  #M4 as Thin Lens
  M4 = SRWLOptL(focLenM4, focLenM4)

  #*****Drift SA -> M4
  SA_M4 = SRWLOptD(pM4)

  #*****Drift M4 -> Sample
  M4_S = SRWLOptD(qM4)

  #*****Propagation Parameters for the Optical Elements
  #Meaning of the array element below:
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

  #Some propagation parameters depend on photon energy
  #        [ 0]  [ 1]  [ 2]  [ 3]  [ 4]  [ 5]  [ 6]  [ 7]  [ 8]  [ 9]  [10]  [11]  [12]  [13]  [14]  [15]  [16]
  pM1_M2_db = {
    250. : [ 0,    0,  1.0,    1,    0,  1.2,  6.0,  1.2,  6.0,    0,    0,    0 ], 
    450. : [ 0,    0,  1.0,    1,    0,  1.2,  4.0,  1.2,  5.0,    0,    0,    0 ],
    1000.: [ 0,    0,  1.0,    1,    0,  1.2,  3.5,  1.2,  3.0,    0,    0,    0 ],
    1500.: [ 0,    0,  1.0,    1,    0,  1.3,  3.0,  1.2,  3.0,    0,    0,    0 ], 
    2000.: [ 0,    0,  1.0,    1,    0,  1.3,  3.0,  1.2,  3.0,    0,    0,    0 ]
  }

  pM2A   = [ 0,    0,  1.0,    0,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pM2_G  = [ 0,    0,  1.0,    1,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pGA    = [ 0,    0,  1.0,    0,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pG     = [ 0,    0,  1.0,    1,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0,   0, sin(defG), cos(defG), 1, 0 ]
  pG_M3  = [ 0,    0,  1.0,    2,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ] #[3]=2 Ensures manipulation with strict values of Rx, Ry
  #pG_M3  = [ 0,    0,  1.0,    1,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ] #[3]=2 Ensures manipulation with strict values of Rx, Ry

  pM3A   = [ 0,    0,  1.0,    0,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pM3    = [ 0,    0,  1.0,    0,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]

  pM3_SA_db = {
    250. : [ 0,    0,  1.0,    4,    0,  1.5,  1.0, 10.0,  1.0,    0,    0,    0 ], #NOTE: reducing resolution may harm...
    450. : [ 0,    0,  1.0,    4,    0,  1.5,  1.0, 10.0,  1.0,    0,    0,    0 ],
    1000.: [ 0,    0,  1.0,    4,    0,  1.0,  1.0,  3.0,  1.0,    0,    0,    0 ],
    1500.: [ 0,    0,  1.0,    4,    0,  1.0,  1.0,  3.0,  1.0,    0,    0,    0 ],
    2000.: [ 0,    0,  1.0,    4,    0,  1.0,  1.0,  3.0,  1.0,    0,    0,    0 ],
  }

  pSA_db = {
    250. : [ 0,    0,  1.0,    0,    0,  0.3,  1.0,  0.1,  2.0,    0,    0,    0 ], 
    450. : [ 0,    0,  1.0,    0,    0,  0.3,  1.0,  0.1,  2.0,    0,    0,    0 ], 
    1000.: [ 0,    0,  1.0,    0,    0,  0.4,  1.0,  0.1,  1.0,    0,    0,    0 ],
    1500.: [ 0,    0,  1.0,    0,    0,  0.4,  1.0,  0.1,  1.0,    0,    0,    0 ], 
    2000.: [ 0,    0,  1.0,    0,    0,  0.4,  1.0,  0.1,  1.0,    0,    0,    0 ]
  }

  pSA_M4 = [ 0,    0,  1.0,    3,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pM4A   = [ 0,    0,  1.0,    0,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pM4    = [ 0,    0,  1.0,    1,    0,  1.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pM4_S  = [ 0,    0,  1.0,    4,    0,  2.0,  1.0,  1.0,  1.0,    0,    0,    0 ]
  pSfin  = [ 0,    0,  1.0,    1,    1,  1.0,  1.5,  0.3,  2.0,    0,    0,    0 ]

  OE  = [  M1_M2,          M2A,  M2_G,  GA,  G,  G_M3,  M3A,  M3,  M3_SA,          SA,          SA_M4,  M4A,  M4,  M4_S]           
  pOE = [ pM1_M2_db[_en], pM2A, pM2_G, pGA, pG, pG_M3, pM3A, pM3, pM3_SA_db[_en], pSA_db[_en], pSA_M4, pM4A, pM4, pM4_S, pSfin]

  #Creating BL "container" (possibly with reduced number of optical elements)
  if((_n_opt_el < 0) or (_n_opt_el >= len(pOE))): return SRWLOptC(OE, pOE)
  elif(_n_opt_el == 0): return None
  else:
    OEr = []; pOEr = []
    for i in range(_n_opt_el):
      OEr.append(OE[i]); pOEr.append(pOE[i])
    return SRWLOptC(OEr, pOEr)

#*********************************Single-Electron Emission and Wavefront Propagation calculation
def calcSE(_e_beam, _mag, _mesh, _bl, _sr_meth=1, _sr_prec=0.01, _sr_samp_fact=1, _fnsuf=''):

  #Wavefront (placeholder)
  wfr = SRWLWfr() 
  wfr.allocate(_mesh.ne, _mesh.nx, _mesh.ny) #Numbers of points vs Photon Energy, Horizontal and Vertical Positions
  wfr.mesh     =  deepcopy(_mesh)
  wfr.partBeam =  deepcopy(_e_beam)

  #Precision Parameters for SR calculation
  #_sr_meth: #SR calculation method: 0- "manual", 1- "auto-undulator", 2- "auto-wiggler"
  #_sr_prec: #relative precision
  zStartInteg = 0. #longitudinal position to start  integration (effective if < zEndInteg)
  zEndInteg   = 0. #longitudinal position to finish integration (effective if > zStartInteg)
  npTraj = 50000 #Number of points for trajectory calculation 
  useTermin = 1 #Use "terminating terms" (i.e. asymptotic expansions at zStartInteg and zEndInteg) or not (1 or 0 respectively)
  #_sr_samp_fact: #sampling factor for adjusting nx, ny (effective if > 0)
  arPrecPar = [_sr_meth, _sr_prec, zStartInteg, zEndInteg, npTraj, useTermin, _sr_samp_fact]

  sys.stdout.write('   Performing Initial Single-E Electric Field calculation ... '); sys.stdout.flush()
  srwl.CalcElecFieldSR(wfr, 0, _mag, arPrecPar)
  sys.stdout.write('done\n')
  mesh0 = deepcopy(wfr.mesh)  
  sys.stdout.write('   Extracting Intensity from the Calculated Initial Electric Field ... '); sys.stdout.flush()
  arI = array('f', [0]*mesh0.nx*mesh0.ny) #"flat" array to take 2D intensity data
  srwl.CalcIntFromElecField(arI, wfr, 6, 0, 3, mesh0.eStart, 0, 0)
  arIx = array('f', [0]*mesh0.nx) #"flat" array to take 1D intensity data
  srwl.CalcIntFromElecField(arIx, wfr, 6, 0, 1, mesh0.eStart, 0, 0)
  arIy = array('f', [0]*mesh0.ny) #"flat" array to take 1D intensity data
  srwl.CalcIntFromElecField(arIy, wfr, 6, 0, 2, mesh0.eStart, 0, 0)
  sys.stdout.write('done\n')
  sys.stdout.write('   Saving the Initial Wavefront Intensity into a file ... '); sys.stdout.flush()
  srwl_uti_save_intens_ascii(arI, mesh0, strIntSE_OutFilePath+_fnsuf)
  sys.stdout.write('done\n')

  sys.stdout.write('   Simulating Electric Field Wavefront Propagation ... '); sys.stdout.flush()
  t0 = time.time()
  srwl.PropagElecField(wfr, _bl)
  sys.stdout.write('done\n   lasted '+repr(round(time.time() - t0))+' s\n')
  mesh1 = deepcopy(wfr.mesh)  
  sys.stdout.write('   Extracting Intensity from the Propagated Electric Field  ... '); sys.stdout.flush()
  arI1 = array('f', [0]*mesh1.nx*mesh1.ny) #"flat" 2D array to take intensity data
  srwl.CalcIntFromElecField(arI1, wfr, 6, 0, 3, mesh1.eStart, 0, 0)
  #srwl.CalcIntFromElecField(arI1, wfr, 0, 5, 3, mesh1.eStart, 0, 0) #ReEx
  #srwl.CalcIntFromElecField(arI1, wfr, 0, 6, 3, mesh1.eStart, 0, 0) #ImEx
  arI1x = array('f', [0]*mesh1.nx) #"flat" array to take 1D intensity data
  srwl.CalcIntFromElecField(arI1x, wfr, 6, 0, 1, mesh1.eStart, 0, 0)
  arI1y = array('f', [0]*mesh1.ny) #"flat" array to take 1D intensity data
  srwl.CalcIntFromElecField(arI1y, wfr, 6, 0, 2, mesh1.eStart, 0, 0)
  sys.stdout.write('done\n')
  sys.stdout.write('   Saving the Propagated Wavefront Intensity data to a file ... '); sys.stdout.flush()
  srwl_uti_save_intens_ascii(arI1, mesh1, strIntPropSE_OutFilePath+_fnsuf)
  sys.stdout.write('done\n')

  sys.stdout.write('   Plotting the results (blocks script execution; close any graph windows to proceed) ... '); sys.stdout.flush()
  plotMesh0x = [1000*mesh0.xStart, 1000*mesh0.xFin, mesh0.nx]
  plotMesh0y = [1000*mesh0.yStart, 1000*mesh0.yFin, mesh0.ny]
  uti_plot2d(arI, plotMesh0x, plotMesh0y, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity Before Propagation'])
  uti_plot1d(arIx, plotMesh0x, ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity (horizontal cut at y = 0)'])
  uti_plot1d(arIy, plotMesh0y, ['Vertical Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity (vertical cut at x = 0)'])
  plotMesh1x = [1000*mesh1.xStart, 1000*mesh1.xFin, mesh1.nx]
  plotMesh1y = [1000*mesh1.yStart, 1000*mesh1.yFin, mesh1.ny]
  uti_plot2d(arI1, plotMesh1x, plotMesh1y, ['Horizontal Position [mm]', 'Vertical Position [mm]', 'Intensity After Propagation'])
  uti_plot1d(arI1x, plotMesh1x, ['Horizontal Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity After Prop. (horizontal cut at y = 0)'])
  uti_plot1d(arI1y, plotMesh1y, ['Vertical Position [mm]', 'Intensity [ph/s/.1%bw/mm^2]', 'Intensity After Prop. (vertical cut at x = 0)'])
  uti_plot_show() #show all graphs (blocks script execution; close all graph windows to proceed)
  sys.stdout.write('all done\n')

#*********************************Multi-Electron Emission and Wavefront Propagation calculation
def calcME(_e_beam, _mag, _mesh, _bl, _bw=0, _sr_meth=1, _sr_prec=0.01, _sr_samp_fact=1, _n_elec=1000000, _n_avg=10, _n_save_per=10, _fnsuf=''):

  ePhInteg = 0 #Switch to integrate over photon energy or not

  if(_bw > 0):
    ePhInteg = 1
    ePhAvg = 0.5*(_mesh.eStart + _mesh.eFin)
    rEph = _bw*ePhAvg #0.03 #Range for the Photon Energy integration [eV]
    _mesh.eStart = ePhAvg - 0.5*rEph
    _mesh.eFin = ePhAvg + 0.5*rEph
    _mesh.ne = 1

  sys.stdout.write('   Starting simulation of Partially-Coherent Wavefront Propagation (takes a lot of time)... ');sys.stdout.flush()
  radStokesProp = srwl_wfr_emit_prop_multi_e(_e_beam, _mag, _mesh, _sr_meth, _sr_prec, _n_elec, _n_avg, _n_save_per, 
                                             strIntPropME_OutFilePath+_fnsuf, _sr_samp_fact, _bl, _e_ph_integ=ePhInteg)
  sys.stdout.write('all done\n')

#*********************************Entry point
if __name__=="__main__":
  import optparse
  p = optparse.OptionParser()
  p.add_option('-i',  '--isamp',  dest='isamp',  metavar="NUMBER",  default=0.07,  type="float",  help="initial wavefront sampling factor")
  p.add_option('-e',  '--energy',  dest='en',  metavar="NUMBER",  default=1000,  type="float",  help="energy [eV]")
  p.add_option('-d',  '--delta',  dest='de',  metavar="NUMBER",  default=0,  type="float",  help="delta energy [eV]")
  p.add_option('-b',  '--bdetune',  dest='bd',  metavar="NUMBER",  default=0,  type="float",  help="relative detuning of undulator magnetic field from resonant value")
  p.add_option('-m',  '--melec',  dest='melec',  metavar="NUMBER",  default=1,  type="int",  help="number of macro-electrons to take into account")
  p.add_option('-a',  '--melavg',  dest='melavg',  metavar="NUMBER",  default=10,  type="int",  help="number of mmacro-electrons (wavefronts) to average on each node (for MPI calculations)")
  p.add_option('-s',  '--melsaveper',  dest='melsaveper',  metavar="NUMBER",  default=10,  type="int",  help="Saving periodicity (in terms of macro-electrons) of the resulting intensity")
  p.add_option('-w',  '--bw',  dest='bw',  metavar="NUMBER",  default=0,  type="float",  help="relative bandwidth (is taken into account at multi-electron calculations)")
  p.add_option('-n',  '--noptel',  dest='noptel',  metavar="NUMBER",  default=-1,  type="int",  help="number of optical elements to propagate wavefront through")
  p.add_option('-f',  '--fnsuf', dest='fnsuf', metavar='FILE', default='', help='output file name suffix')
  opt, args = p.parse_args()

  print(help_str)
  mag = setUndulator(opt.en, opt.bd)
  enPhExact = opt.en + opt.de
  mesh = setRadMesh(enPhExact)
  e_beam = setElecBeam(-1.9)
  beamline = setBeamline(opt.en, opt.noptel)

  fileNameSuf = opt.fnsuf+'_'+repr(enPhExact)+"eV.dat"

  if(opt.melec <= 1):
    calcSE(e_beam, mag, mesh, beamline, _sr_samp_fact=opt.isamp, _fnsuf=fileNameSuf)
  else:
    calcME(e_beam, mag, mesh, beamline, _bw=opt.bw, _sr_samp_fact=opt.isamp,
           _n_elec=opt.melec, _n_avg=opt.melavg, _n_save_per=opt.melsaveper, _fnsuf=fileNameSuf)

