#############################################################################
# uti_math module: misc. mathematical utilities / functions
# v 0.01
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from array import *
from math import *
from copy import *
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
#def interp_1d_lin_var(_x, _ar_x, _ar_f):
def interp_1d_var(_x, _ar_x, _ar_f, _ord=3):
    """
    Interpolate linearly 1D function value tabulated on non-equidistant (irregular) mesh
    :param _x: argument at which function value should be calculated
    :param _ar_x: array or list of increasing argument values, at which the function is tabulated
    :param _ar_f: array or list of tabulated values corresponding to arguments in _ar_x
    :param _ord: order of polynomial interpolation (1- linear, 2- quadratic, 3- cubic)    
    :return: function value found by polynomial interpolation
    """

    sErrBadArrays = 'Incorrect/incompatible lengths of argument and function value arrays.'

    nx = len(_ar_x)
    if(nx <= 0): raise Exception(sErrBadArrays)
    nf = len(_ar_f)
    if(nf <= 0): raise Exception(sErrBadArrays)
    
    if(nx > nf): nx = nf
    if(nx < 2): raise Exception(sErrBadArrays)
    
    nx_mi_1 = nx - 1

    if(_x <= _ar_x[0]): return _ar_f[0]
    elif(_x >= _ar_x[nx_mi_1]): return _ar_f[nx_mi_1]

    if(_ord > nx_mi_1): _ord = nx_mi_1

    i0 = 0
    for i in range(1, nx):
        if(_x < _ar_x[i]):
            i0 = i - 1
            break

    if(_ord == 1):
        #return interp_1d(_x, _ar_x[i0], _ar_x[i0+1] - _ar_x[i0], 2, [_ar_f[i0], _ar_f[i0+1]], 1)
        x0 = _ar_x[i0]
        step = _ar_x[i0+1] - x0
        t = (_x - x0)/step
        f0 = _ar_f[i0]
        return f0 + (_ar_f[i0+1] - f0)*t

    elif(_ord == 2):
        im1 = i0 - 1
        ip1 = i0 + 1
        #nm1 = nx - 1
        if(ip1 < 0):
            im1 = 0
            i0 = 1
            ip1 = 2
        elif(ip1 > nx_mi_1):
            im1 = nx - 3
            i0 = nx - 2
            ip1 = nx_mi_1

        xm1 = _ar_x[im1]
        x0 = _ar_x[i0]
        xp1 = _ar_x[ip1]
        dxm1 = abs(_x - xm1)
        dx0 = abs(_x - x0)
        dxp1 = abs(_x - xp1)
        if(dxm1 < dx0):
            if(im1 > 0):
                im1 -= 1
                i0 -= 1
                ip1 -= 1
        elif(dxp1 < dx0):
            if(ip1 < nx_mi_1):
                im1 += 1
                i0 += 1
                ip1 += 1

        x0 = _ar_x[i0]
        dxm1 = _ar_x[im1] - x0
        dxp1 = _ar_x[ip1] - x0
        fm1 = _ar_f[im1]
        f0 = _ar_f[i0]
        fp1 = _ar_f[ip1]

        invD = 1./(dxm1* dxp1*(dxm1 - dxp1))
        a = ((dxm1 - dxp1)*f0 + dxp1*fm1 - dxm1*fp1)*invD
        b = (dxp1*dxp1*(f0 - fm1) + dxm1*dxm1*(fp1 - f0))*invD
        dx = _x - x0
        return (a*dx + b)*dx + f0

    elif(_ord == 3):
        im1 = i0 - 1
        ip1 = i0 + 1
        ip2 = i0 + 2
        if(im1 < 0):
            im1 = 0
            i0 = 1
            ip1 = 2
            ip2 = 3
        elif(ip2 > nx_mi_1):
            im1 = nx - 4
            i0 = nx - 3
            ip1 = nx - 2
            ip2 = nx_mi_1
            
        x0 = _ar_x[i0]
        dxm1 = _ar_x[im1] - x0
        dxp1 = _ar_x[ip1] - x0
        dxp2 = _ar_x[ip2] - x0
        fm1 = _ar_f[im1]
        f0 = _ar_f[i0]
        fp1 = _ar_f[ip1]
        fp2 = _ar_f[ip2]
        #print(_x - x0, dxm1, dxp1, dxp2)
        #print(fm1, f0, fp1, fp2)

        invD = 1./(dxm1*dxp1*dxp2*(dxm1 - dxp1)*(dxm1 - dxp2)*(dxp1 - dxp2))
        invD1 = 1./(dxm1*dxp1*dxp2)
        dxm1e2 = dxm1*dxm1
        dxm1e3 = dxm1e2*dxm1
        dxp1e2 = dxp1*dxp1
        dxp1e3 = dxp1e2*dxp1
        dxp2e2 = dxp2*dxp2
        dxp2e3 = dxp2e2*dxp2
        a1 = (-dxp1e2*(dxp1 - dxp2)*dxp2e2*(f0 - fm1) + dxm1e2*(dxp2e3*(fp1 - f0) + dxp1e3*(f0 - fp2)) + dxm1e3*(dxp2e2*(f0 - fp1) + dxp1e2*(fp2 - f0)))*invD
        a2 = ((dxm1 + dxp1 + dxp2)*f0)*invD1 + (dxm1*dxp2*(dxm1e2 - dxp2e2)*fp1 + dxp1e3*(dxm1*fp2 - dxp2*fm1) + dxp1*(dxp2e3*fm1 - dxm1e3*fp2))*invD
        a3 = -f0*invD1 + (dxm1*dxp2*(dxp2 - dxm1)*fp1 + dxp1e2*(dxp2*fm1 - dxm1*fp2) + dxp1*(dxm1e2*fp2 - dxp2e2*fm1))*invD
        dx = _x - x0
        return ((a3*dx + a2)*dx + a1)*dx + f0

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
##    if(_x == 0.): return _x
##    sgn = 1.
##    if(_x < 0.):
##        _x = -_x
##        sgn = -1
##    order = round(log10(_x))
##    fact = 10**order
##    roundNum = round(_x/fact, _ndig)
##    res = roundNum*fact*sgn
    res = round(_x, _ndig)
    return res

#****************************************************************************
def integ_array(_ar, _h, _dupl=False):
    """
    Integrates array (or list), eventually making a copy of it before the integration
    :param _ar: array to integrate
    :param _h: step size
    :param _dupl: duplicate the magnetic field object or not
    """
    ar = None
    if(_dupl): ar = deepcopy(_ar)
    else: ar = _ar

    hd2 = 0.5*_h
    auxInt = 0
    lenAr = len(_ar)
    for i in range(lenAr - 1):
        ar_i = _ar[i] #in case if array has to be integrated in place
        ar[i] = auxInt
        auxInt += hd2*(ar_i + _ar[i + 1])
    ar[lenAr - 1] = auxInt
    return ar

#****************************************************************************
def matr_prod(_A, _B):
    """
    Multiplies matrix _A by matrix _B 
    """
    # Matrix multiplication
    B0 = _B[0]
    lenB = len(_B)
    lenA = len(_A)
    if(len(_A[0]) != lenB): # Check matrix dimensions        
        Exception('Matrices have wrong dimensions')
    if(isinstance(B0, list) or isinstance(B0, array) or isinstance(B0, tuple)): #_B is matrix
        lenB0 = len(B0)
        C = [[0 for row in range(lenB0)] for col in range(lenA)]
        for i in range(lenA):
            for j in range(lenB0):
                for k in range(lenB):
                    C[i][j] += _A[i][k]*_B[k][j]
    else: #_B is vector
        C = [0 for row in range(lenB)]
        for i in range(lenA):
            for k in range(lenB):
                C[i] += _A[i][k]*_B[k]
    return C

#****************************************************************************
def matr_print(_A):
    """
    Prints matrix _A
    """
    for i in range(len(_A)):
        print(_A[i])
 
#****************************************************************************
def trf_rotation(_V, _ang, _P):
    """
    Sets up matrix and vector describing rotation about axis _V passing through a point _P about an angle _ang
    :param _V: vector (array of 3 Cartesian coordinates) rdefining rotation axis
    :param _ang: rotation angle [rad]
    :param _P: point (array of 3 Cartesian coordinates) rotation axis passes through
    :returns list containing the 3x3 matrix and 3-element vector
    """
    normFact = 1./sqrt(_V[0]*_V[0] + _V[1]*_V[1] + _V[2]*_V[2]);
    axVect = [normFact*_V[0], normFact*_V[1], normFact*_V[2]]
    VxVx = axVect[0]*axVect[0]
    VyVy = axVect[1]*axVect[1]
    VzVz = axVect[2]*axVect[2]
    cosAng = cos(_ang)
    sinAng = sin(_ang)
    one_m_cos = 1. - cosAng
    one_m_cosVxVy = one_m_cos*axVect[0]*axVect[1]
    one_m_cosVxVz = one_m_cos*axVect[0]*axVect[2]
    one_m_cosVyVz = one_m_cos*axVect[1]*axVect[2]
    sinVx = sinAng*axVect[0]
    sinVy = sinAng*axVect[1]
    sinVz = sinAng*axVect[2]
    st0 = [VxVx + cosAng*(VyVy + VzVz), one_m_cosVxVy - sinVz, one_m_cosVxVz + sinVy]
    st1 = [one_m_cosVxVy + sinVz, VyVy + cosAng*(VxVx + VzVz), one_m_cosVyVz - sinVx]
    st2 = [one_m_cosVxVz - sinVy, one_m_cosVyVz + sinVx, VzVz + cosAng*(VxVx + VyVy)]
    M = [st0, st1, st2]
    st00 = [1. - st0[0], -st0[1], -st0[2]]
    st01 = [-st1[0], 1. - st1[1], -st1[2]]
    st02 = [-st2[0], -st2[0], 1. - st2[2]]
    M0 = [st00, st01, st02]
    V = matr_prod(M0, _P)
    return [M, V]
