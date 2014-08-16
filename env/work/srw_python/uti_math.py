#############################################################################
# uti_math module: misc. mathematical utilities / functions
# v 0.01
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from array import *
from math import *
#from copy import *
#import random
#import sys
#import os

#****************************************************************************
def interp_1d(_x, _x_min, _x_step, _nx, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
    """
    Interpolate 1D function value tabulated on equidistant mesh, using polynomial interpolation
    :param _x: argument at which function value should be calculated
    :param _x_min: minimal argument value of the tabulated function
    :param _x_step: step of mesh at which function is tabulated
    :param _nx: number of points in mesh at which function is tabulated
    :param _ar_f: tabulated function list or array
    :param _ord: order of polynomial interpolation (1- linear, 2- quadratic, 3- cubic)
    :param _ix_per: argument index period of function data alignment (e.g. to interpolate one component of complex data, or in one dimension of multi-dimensional data)
    :param _ix_ofst: argument index offset of function data alignment
    :return: function value found by polynomial interpolation
    """
    if(_ord == 1):
        i0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(i0 < 0):
            i0 = 0
        elif(i0 >= _nx - 1):
            i0 = _nx - 2
        i1 = i0 + 1
        f0 = _ar_f[i0*_ix_per + _ix_ofst]
        f1 = _ar_f[i1*_ix_per + _ix_ofst]
        t = (_x - (_x_min + _x_step*i0))/_x_step
        return f0 + (f1 - f0)*t
    elif(_ord == 2):
        i0 = int(round((_x - _x_min)/_x_step))
        if(i0 < 1):
            i0 = 1
        elif(i0 >= _nx - 1):
            i0 = _nx - 2
        im1 = i0 - 1
        i1 = i0 + 1
        t = (_x - (_x_min + _x_step*i0))/_x_step
        a0 = _ar_f[i0*_ix_per + _ix_ofst]
        fm1 = _ar_f[im1*_ix_per + _ix_ofst]
        f1 = _ar_f[i1*_ix_per + _ix_ofst]
        a1 = 0.5*(f1 - fm1)
        a2 = 0.5*(fm1 + f1 - 2*a0)
        return a0 + t*(a1 + t*a2)
    elif(_ord == 3):
        i0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(i0 < 1):
            i0 = 1
        elif(i0 >= _nx - 2):
            i0 = _nx - 3
        im1 = i0 - 1
        i1 = i0 + 1
        i2 = i0 + 2
        t = (_x - (_x_min + _x_step*i0))/_x_step
        a0 = _ar_f[i0*_ix_per + _ix_ofst]
        fm1 = _ar_f[im1*_ix_per + _ix_ofst]
        f1 = _ar_f[i1*_ix_per + _ix_ofst]
        f2 = _ar_f[i2*_ix_per + _ix_ofst]
        a1 = -0.5*a0 + f1 - f2/6. - fm1/3.
        a2 = -a0 + 0.5*(f1 + fm1)
        a3 = 0.5*(a0 - f1) + (f2 - fm1)/6.
        return a0 + t*(a1 + t*(a2 + t*a3))
    return 0

#****************************************************************************
def interp_2d(_x, _y, _x_min, _x_step, _nx, _y_min, _y_step, _ny, _ar_f, _ord=3, _ix_per=1, _ix_ofst=0):
    """
    Interpolate 2D function value tabulated on equidistant rectangular mesh and represented by C-aligned flat array, using polynomial interpolation
    :param _x: first argument at which function value should be calculated
    :param _y: second argument at which function value should be calculated
    :param _x_min: minimal value of the first argument of the tabulated function
    :param _x_step: step of the first argument at which function is tabulated
    :param _nx: number of points vs first argument at which function is tabulated
    :param _y_min: minimal value of the second argument of the tabulated function
    :param _y_step: step of the second argument at which function is tabulated
    :param _ny: number of points vs second argument at which function is tabulated
    :param _ar_f: function tabulated on 2D mesh, aligned as "flat" C-type list or array (first argument is changing most frequently)
    :param _ord: "order" of polynomial interpolation (1- bi-linear (on 4 points), 2- "bi-quadratic" (on 6 points), 3- "bi-cubic" (on 12 points))
    :param _ix_per: period of first argument index of the function data alignment (e.g. to interpolate one component of complex data, or in one dimension of multi-dimensional data)
    :param _ix_ofst: offset of the first argument index in function data alignment
    :return: function value found by 2D polynomial interpolation
    """
    if(_ord == 1): #bi-linear interpolation based on 4 points
        ix0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(ix0 < 0):
            ix0 = 0
        elif(ix0 >= _nx - 1):
            ix0 = _nx - 2
        ix1 = ix0 + 1
        tx = (_x - (_x_min + _x_step*ix0))/_x_step
        
        iy0 = int(trunc((_y - _y_min)/_y_step + 1.e-09))
        if(iy0 < 0):
            iy0 = 0
        elif(iy0 >= _ny - 1):
            iy0 = _ny - 2
        iy1 = iy0 + 1
        ty = (_y - (_y_min + _y_step*iy0))/_y_step

        nx_ix_per = _nx*_ix_per
        iy0_nx_ix_per = iy0*nx_ix_per
        iy1_nx_ix_per = iy1*nx_ix_per
        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        a10 = f10 - a00
        a01 = f01 - a00
        a11 = a00 - f01 - f10 + f11
        return a00 + tx*(a10 + ty*a11) + ty*a01

    elif(_ord == 2): #bi-quadratic interpolation based on 6 points
        ix0 = int(round((_x - _x_min)/_x_step))
        if(ix0 < 1):
            ix0 = 1
        elif(ix0 >= _nx - 1):
            ix0 = _nx - 2
        ixm1 = ix0 - 1
        ix1 = ix0 + 1
        tx = (_x - (_x_min + _x_step*ix0))/_x_step

        iy0 = int(round((_y - _y_min)/_y_step))
        if(iy0 < 1):
            iy0 = 1
        elif(iy0 >= _ny - 1):
            iy0 = _ny - 2
        iym1 = iy0 - 1
        iy1 = iy0 + 1
        ty = (_y - (_y_min + _y_step*iy0))/_y_step

        nx_ix_per = _nx*_ix_per
        iym1_nx_ix_per = iym1*nx_ix_per
        iy0_nx_ix_per = iy0*nx_ix_per
        iy1_nx_ix_per = iy1*nx_ix_per
        ixm1_ix_per_p_ix_ofst = ixm1*_ix_per + _ix_ofst
        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
        fm10 = _ar_f[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f0m1 = _ar_f[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        a10 = 0.5*(f10 - fm10)
        a01 = 0.5*(f01 - f0m1)
        a11 = a00 - f01 - f10 + f11
        a20 = 0.5*(f10 + fm10) - a00
        a02 = 0.5*(f01 + f0m1) - a00
        return a00 + tx*(a10 + tx*a20 + ty*a11) + ty*(a01 + ty*a02)
    
    elif(_ord == 3): #bi-cubic interpolation based on 12 points
        ix0 = int(trunc((_x - _x_min)/_x_step + 1.e-09))
        if(ix0 < 1):
            ix0 = 1
        elif(ix0 >= _nx - 2):
            ix0 = _nx - 3
        ixm1 = ix0 - 1
        ix1 = ix0 + 1
        ix2 = ix0 + 2
        tx = (_x - (_x_min + _x_step*ix0))/_x_step

        iy0 = int(trunc((_y - _y_min)/_y_step + 1.e-09))
        if(iy0 < 1):
            iy0 = 1
        elif(iy0 >= _ny - 2):
            iy0 = _ny - 3
        iym1 = iy0 - 1
        iy1 = iy0 + 1
        iy2 = iy0 + 2
        ty = (_y - (_y_min + _y_step*iy0))/_y_step

        nx_ix_per = _nx*_ix_per
        iym1_nx_ix_per = iym1*nx_ix_per
        iy0_nx_ix_per = iy0*nx_ix_per
        iy1_nx_ix_per = iy1*nx_ix_per
        iy2_nx_ix_per = iy2*nx_ix_per
        ixm1_ix_per_p_ix_ofst = ixm1*_ix_per + _ix_ofst
        ix0_ix_per_p_ix_ofst = ix0*_ix_per + _ix_ofst
        ix1_ix_per_p_ix_ofst = ix1*_ix_per + _ix_ofst
        ix2_ix_per_p_ix_ofst = ix2*_ix_per + _ix_ofst
        f0m1 = _ar_f[iym1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f1m1 = _ar_f[iym1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        fm10 = _ar_f[iy0_nx_ix_per + ixm1_ix_per_p_ix_ofst]
        a00 = _ar_f[iy0_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f10 = _ar_f[iy0_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f20 = _ar_f[iy0_nx_ix_per + ix2_ix_per_p_ix_ofst]
        fm11 = _ar_f[iy1_nx_ix_per + ixm1_ix_per_p_ix_ofst]
        f01 = _ar_f[iy1_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f11 = _ar_f[iy1_nx_ix_per + ix1_ix_per_p_ix_ofst]
        f21 = _ar_f[iy1_nx_ix_per + ix2_ix_per_p_ix_ofst]
        f02 = _ar_f[iy2_nx_ix_per + ix0_ix_per_p_ix_ofst]
        f12 = _ar_f[iy2_nx_ix_per + ix1_ix_per_p_ix_ofst]
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
        return a00 + tx*(a10 + tx*(a20 + tx*(a30 + ty*a31) + ty*a21) + ty*a11) + ty*(a01 + ty*(a02 + ty*(a03 + tx*a13) + tx*a12))
    return 0

#****************************************************************************
def num_round(_x, _ndig=8):
    if(_x == 0.): return _x
    sgn = 1
    if(_x < 0.):
        _x = -_x
        sgn = -1
    order = round(log10(_x))
    fact = 10**order
    return round(_x/fact, _ndig)*fact*sgn
