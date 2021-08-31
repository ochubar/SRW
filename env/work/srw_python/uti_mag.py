
#############################################################################
# Magnetic Field Utilities Module (to be used e.g. with Radia code)
# v 0.01
#############################################################################

from __future__ import absolute_import, division, print_function #Py 2.*/3.* compatibility
from math import *
from array import *

#***************Harmonic Decomposition of Periodic Magnetic Field
def uti_mag_fld_harm(_fld_over_per, _n=1):

    nFld = len(_fld_over_per)
    nn = nFld*_n #21*_n
    nn_mi_1 = nn - 1

    angStep = 2*pi*_n/nn_mi_1
    ang = 0.
    sumC = 0.; sumS = 0.
    
    for i in range(nn):
        cosAng = cos(ang); sinAng = sin(ang)
        fc = cosAng*_fld_over_per[i]
        fs = sinAng*_fld_over_per[i]
        if((i == 0) or (i == nn_mi_1)):
            fc *= 0.5; fs *= 0.5
        sumC += fc; sumS += fs
        ang += angStep

    mult = 2./nn_mi_1
    sumC *= mult; sumS *= mult

    #OC15042021
    ang = 0
    if(sumC > 0):
        ang = atan(sumS/sumC)
    elif(sumC == 0):
        if(sumS >= 0): ang = 0.5*pi
        else: ang = -0.5*pi
    else:
        if(sumS >= 0):
            ang = pi - atan(abs(sumS/sumC))
        else:
            ang = pi + atan(abs(sumS/sumC))
    return [sqrt(sumC*sumC + sumS*sumS), ang]

    #if(sumC != 0): return [sqrt(sumC*sumC + sumS*sumS), ang]
    #else: return [sqrt(sumC*sumC + sumS*sumS), pi/2] #OC08042021 (consider different cases!)
    #return [sqrt(sumC*sumC + sumS*sumS), atan(sumS/sumC)]

#***************Deflection Parameter of Undulator / Wiggler
def uti_und_keff(_fld_over_per, _per_mm):

    nDim = 1
    len_fld_over_per = len(_fld_over_per)
    if((len_fld_over_per == 2) and (isinstance(_fld_over_per[0], list) or isinstance(_fld_over_per[0], array))): nDim = 2

    Keff = 0.
    if(nDim == 1):
        B0 = uti_mag_fld_harm(_fld_over_per)[0]
        Keff = 0.0933729*_per_mm*B0
    elif(nDim == 2):
        B01 = uti_mag_fld_harm(_fld_over_per[0])[0]
        B02 = uti_mag_fld_harm(_fld_over_per[1])[0]
        Keff = 0.0933729*_per_mm*sqrt(B01*B01 + B02*B02)
        
    return Keff

#***************Fundamental Photon Energy of Undulator Radiation
def uti_und_e1_from_fld(_fld_over_per, _per_mm, _el_en_GeV):
    Keff = uti_und_keff(_fld_over_per, _per_mm)
    return 9496.3422*_el_en_GeV*_el_en_GeV/_per_mm/(1. + 0.5*Keff*Keff)

def uti_und_e1_from_k(_k, _per_mm, _el_en_GeV):
    return 9496.3422*_el_en_GeV*_el_en_GeV/_per_mm/(1. + 0.5*_k*_k)
   

