#############################################################################
# SRWLib for Python Utility module for Synchrotron Sources (electron beams and ID parameters)
# Author: O.C.
# v 0.01
#############################################################################

from srwlib import *

#****************************************************************************
def srwl_uti_src_e_beam(_nm):
    """Instantiates electron beam structures describing different existing sources
    :param _name: string identifying a source
    :return: SRWLPartBeam object
    """

    #E-Beam params in the order: _Iavg, _e, _sig_e, _emit_x, _beta_x, _alpha_x, _eta_x, _eta_x_pr, _emit_y, _beta_y, _alpha_y
    allBeams = [
        ['NSLS-II Low Beta Day 1',  [0.5,    3,  0.89e-03,   0.9e-09,  2.02,     0,     0,      0,     8e-12,   1.06,      0]],
        ['NSLS-II Low Beta Final',  [0.5,    3,  0.89e-03,  0.55e-09,  2.02,     0,     0,      0,     8e-12,   1.06,      0]],
        ['NSLS-II High Beta Day 1', [0.5,    3,  0.89e-03,   0.9e-09, 20.85,     0,     0,      0,     8e-12,    3.4,      0]],
        ['NSLS-II High Beta Final', [0.5,    3,  0.89e-03,  0.55e-09, 20.85,     0,     0,      0,     8e-12,    3.4,      0]],
        ['NSLS-II 3PW Day 1',       [0.5,    3,  0.89e-03,   0.9e-09, 2.956, 1.932, 0.137, -0.105,     8e-12, 19.653, -0.806]],
        ['NSLS-II 3PW Final',       [0.5,    3,  0.89e-03,  0.55e-09, 2.956, 1.932, 0.137, -0.105,     8e-12, 19.653, -0.806]],
        ['NSLS-II BM Day 1',        [0.5,    3,  0.89e-03,   0.9e-09,   1.5,     0, 0.137,   -0.1,     8e-12,   22.5,   -0.9]],
        ['NSLS-II BM Final',        [0.5,    3,  0.89e-03,  0.55e-09,   1.5,     0, 0.137,   -0.1,     8e-12,   22.5,   -0.9]],
        ['DIAMOND Low Beta',        [0.3,    3,   1.0e-03, 2.797e-09, 2.283,     0,     0,      0, 24.32e-12,  2.514,      0]],
        ['SOLEIL Short',            [0.5, 2.75, 1.016e-03,  3.73e-09, 17.78,     0,  0.28,      0,    37e-12,   1.75,      0]],
        ['SOLEIL Medium',           [0.5, 2.75, 1.016e-03,  3.73e-09,   4.0,     0,  0.13,      0,    37e-12,   1.77,      0]],
        ['SOLEIL Long',             [0.5, 2.75, 1.016e-03,  3.73e-09, 10.09,     0,   0.2,      0,    37e-12,   8.01,      0]],
        ['SOLEIL BM 1 deg.',        [0.5, 2.75, 1.016e-03,  3.73e-09, 0.603, 0.776, 0.039, -0.088,    37e-12,  16.53,  0.931]],
        ['SOLEIL BM 4 deg.',        [0.5, 2.75, 1.016e-03,  3.73e-09, 0.375, 0.024, 0.021, -0.037,    37e-12,  16.01,  0.899]],
        ['ESRF Low Beta',           [0.2, 6.03,     1e-03,     4e-09,  0.35,     0, 0.031,      0,     5e-12,   2.97,      0]],
        ['ESRF High Beta',          [0.2, 6.03,     1e-03,     4e-09,  37.6,     0, 0.134,      0,     5e-12,   2.95,      0]],
        ['APS',                     [0.1,    7,  0.96e-03,  2.79e-09,  22.7,     0, 0.206,      0,   8.4e-12,    3.1,      0]],
        ['SPring8 High Beta',       [0.1,    8,   1.1e-03,   3.4e-09,  22.6,     0, 0.107,      0,   6.8e-12,    5.6,      0]]
    ]#add more beams

    sTest = _nm.replace(' ', '')
    sTest = sTest.replace('-', '')
    sTest = sTest.capitalize()

    resBeam = SRWLPartBeam()
    nBeams = len(allBeams)
    for i in range(nBeams):
        curInf = allBeams[i]
        curStr = curInf[0]
        curStr = curStr.replace(' ', '')
        curStr = curStr.replace('-', '')
        curStr = curStr.capitalize()
        if sTest == curStr:
            ar = curInf[1]
            resBeam.from_Twiss(_Iavg=ar[0], _e=ar[1], _sig_e=ar[2], _emit_x=ar[3], _beta_x=ar[4], _alpha_x=ar[5], _eta_x=ar[6], _eta_x_pr=ar[7], _emit_y=ar[8], _beta_y=ar[9], _alpha_y=ar[10])
            return resBeam
    return None

