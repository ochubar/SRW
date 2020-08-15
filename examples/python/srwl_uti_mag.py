#############################################################################
# SRWLib for Python: Magnet Utilities v 0.02
#############################################################################

from srwlib import *
#from copy import *

#****************************************************************************
def srwl_mag_kick(_el_en=3., _ang=1., _x_or_y='x', _len=1., _led=0):
    """Setup dipole "kick" magnet
    :param _el_en: electron energy [GeV]
    :param _ang: kick angle [rad]
    :param _x_or_y: plane (horizontal or vertical)
    :param _len: effective magnet length [m]
    :param _led: "soft" edge length for field variation from 10% to 90% [m]; G/(1 + ((z-zc)/d)^2)^2 fringe field dependence is assumed [m]
    :param _zc: longitudinal position of the magnet center
    :param _add: add (=1) or reset (=0) the new magnet structure to the existing approximate magnetic field container
    """

    if(_el_en <= 0): raise Exception('Electron Beam structure is not defined (it is required for defining kick magnet parameters from angle)')
    if(_len <= 0): raise Exception('Inconsistent input magnet parameters: effective length can not be negative')

    if(_led > 0.):
        d = 0.8078259211948791*_led #d parameter in G/(1 + ((z-zc)/d)^2)^2
        L0 = _len - 1.5707963267948966*d #length of const. magnetic field part
        if(L0 < 0.): raise Exception('Inconsistent input magnet parameters: magnet edge length if too large for given effective length')

    B = -3.33564095*_ang*_el_en/_len #sign to be checked!
    n_or_s = 'n'
    if(_x_or_y == 'y'): n_or_s = 's' #to check
    
    return SRWLMagFldM(_G=B, _m=1, _n_or_s=n_or_s, _Leff=_len, _Ledge=_led)

#****************************************************************************
def srwl_mag_extrap_grad_off_mid_plane(_mag_mid, _ry, _ny, _grad_mult=1):
    """
    Extrapolates magnetic field off horizontal mid-plane based on field gradient in the mid-plane (e.g. for dipole with gradient or for a guad).
    :param _mag_mid: input tabulated 3D magnetic field on 2D mesh vs x and z in the horizontal mid-plane in a container (object of SRWLMagFldC type)
    :param _ry: vertical position range of the final 3D mesh [m]
    :param _ny: number of points vs vertical position in the final 3D mesh
    :param _grad_mult: a number the gradient has to be multiplied by
    :returns: resulting exrapolated 3D magnetic field structure in a container (object of SRWLMagFldC type)
    """

    if((_mag_mid is None) or (_mag_mid.arMagFld is None) or (len(_mag_mid.arMagFld) != 1) or (_ry < 0.) or (_ny <= 1)):
        raise Exception("Incorrect input parameters for magnetic field to be extrapolated")

    fld3d_mid = _mag_mid.arMagFld[0]
    #arBx_mid = fld3d_mid.arBx
    arBy_mid = fld3d_mid.arBy
    arBz_mid = fld3d_mid.arBz

    xc = _mag_mid.arXc[0]
    xStart = xc - 0.5*fld3d_mid.rx
    nx = fld3d_mid.nx
    xStep = fld3d_mid.rx/(nx - 1)

    yc = 0
    yStart = -0.5*_ry
    ny = int(_ny)
    yStep = _ry/(ny - 1)

    zc = _mag_mid.arZc[0]
    zStart = zc - 0.5*fld3d_mid.rz
    nz = fld3d_mid.nz
    zStep = fld3d_mid.rz/(nz - 1)

    nTotRes = int(nx*ny*nz)
    arBxRes = array('d', [0]*nTotRes)
    arByRes = array('d', [0]*nTotRes)
    arBzRes = None if(arBz_mid is None) else array('d', [0]*nTotRes)
    
    nx_mi_1 = nx - 1
    #print(nx, xStart, xStep, ny, yStart, yStep, nz, zStart, zStep)
    
    for iz in range(nz):
        iz_nx = iz*nx
        iz_nx_ny = iz*nx*ny
        for ix in range(nx):
            b1y = 0; b2y = 0; b0y = 0; dx = xStep
            if(ix == 0):
                b1y = arBy_mid[iz_nx]
                b2y = arBy_mid[iz_nx + 1]
                b0y = b1y
            elif(ix == nx_mi_1):
                b1y = arBy_mid[iz_nx + nx_mi_1 - 1]
                b2y = arBy_mid[iz_nx + nx_mi_1]
                b0y = b2y
            else:
                b1y = arBy_mid[iz_nx + ix - 1]
                b0y = arBy_mid[iz_nx + ix]
                b2y = arBy_mid[iz_nx + ix + 1]
                dx = 2*xStep
                
            curGrad = _grad_mult*(b2y - b1y)/dx
            y = yStart
            for iy in range(ny):
                ofst = iz_nx_ny + iy*nx + ix
                arBxRes[ofst] = y*curGrad
                arByRes[ofst] = b0y
                y += yStep

    fld3dRes = SRWLMagFld3D(arBxRes, arByRes, arBzRes, nx, ny, nz, fld3d_mid.rx, _ry, fld3d_mid.rz)
    #print(fld3dRes.nx, fld3dRes.rx, fld3dRes.ny, fld3dRes.ry, fld3dRes.nz, fld3dRes.rz)

    return SRWLMagFldC(fld3dRes, xc, yc, zc)

#****************************************************************************
def srwl_mag_track_e_beam_mom(_e_beam, _mag, _z_or_ct, _ct_start=0, _ct_end=0, _npart=1000, _np_in_trj=10000):
    """Estimate Electron Beam Statistical Moments at a given long. pos. (for tests)
    :param _e_beam: input electron beam structure (object of SRWLPartBeam type)
    :param _mag: input magnetic field container (object of SRWLMagFldC type)
    :param _z_or_ct: longitudinal position [m] where the stat. moments have to be calculated, or c*t (if _ct_start != _ct_end)
    :param _ct_start: start time ct moment [m] for trajectory calculation
    :param _ct_end: end timect  moment [m] for trajectory calculation
    :param _npart: number of macro-particles
    :param _np_in_trj: number of points in trajectory
    """

    partTraj = SRWLPrtTrj()
    partTraj.allocate(_np_in_trj, True)
    partTraj.ctStart = _ct_start #Start Time for the calculation
    partTraj.ctEnd = _ct_end #magFldCnt.arMagFld[0].rz

    x0 = _e_beam.partStatMom1.x
    xp0 = _e_beam.partStatMom1.xp
    y0 = _e_beam.partStatMom1.y
    yp0 = _e_beam.partStatMom1.yp
    #en0GeV = _e_beam.partStatMom1.get_E()
    gamma0 = _e_beam.partStatMom1.gamma

    sigXe2 = _e_beam.arStatMom2[0] #[0]: <(x-x0)^2>
    mXXp = _e_beam.arStatMom2[1] #[1]: <(x-x0)*(xp-xp0)>
    sigXpe2 = _e_beam.arStatMom2[2] #[2]: <(xp-xp0)^2>
    sigYe2 = _e_beam.arStatMom2[3] #[3]: <(y-y0)^2>
    mYYp = _e_beam.arStatMom2[4] #[4]: <(y-y0)*(yp-yp0)>
    sigYpe2 = _e_beam.arStatMom2[5] #[5]: <(yp-yp0)^2>
    relEnSpr = sqrt(_e_beam.arStatMom2[10]) #<(E-E0)^2>/E0^2
    #absEnSpr = en0GeV*relEnSpr

    multX = 0.5/(sigXe2*sigXpe2 - mXXp*mXXp)
    BX = sigXe2*multX
    GX = sigXpe2*multX
    AX = mXXp*multX
    sigPX = 1/sqrt(2*GX)
    sigQX = sqrt(sigXpe2)

    multY = 0.5/(sigYe2*sigYpe2 - mYYp*mYYp)
    BY = sigYe2*multY
    GY = sigYpe2*multY
    AY = mYYp*multY
    sigPY = 1/sqrt(2*GY)
    sigQY = sqrt(sigYpe2)

    randAr = array('d', [0]*5) #for random Gaussian numbers
    random.seed(12345)

    ind0 = 0
    r0 = 0
    ctStep = 0
    treatCT = False
    np_in_trj_mi_1 = _np_in_trj - 1
    if(_ct_start != _ct_end): #treat _z_or_ct as ct
        ctStep = (_ct_end - _ct_start)/np_in_trj_mi_1
        ind0 = int((_z_or_ct - _ct_start)/ctStep)
        if((ind0 >= 0) and (ind0 < np_in_trj_mi_1)): 
            r0 = (_z_or_ct - (_ct_start + ind0*ctStep))/ctStep
            treatCT = True
        if(treatCT is not True):
            raise Exception("Trajectory point corresponding to given time moment was not found")

    xSum = 0; ySum = 0
    xpSum = 0; ypSum = 0
    xe2Sum = 0; ye2Sum = 0
    xpe2Sum = 0; ype2Sum = 0
    xxpSum = 0; yypSum = 0
    for i in range(_npart):

        for ir in range(5):
            randAr[ir] = random.gauss(0, 1)

        auxPXp = sigQX*randAr[0]
        auxPX = sigPX*randAr[1] + AX*auxPXp/GX
        partTraj.partInitCond.x = x0 + auxPX
        partTraj.partInitCond.xp = xp0 + auxPXp
        auxPYp = sigQY*randAr[2]
        auxPY = sigPY*randAr[3] + AY*auxPYp/GY
        partTraj.partInitCond.y = y0 + auxPY
        partTraj.partInitCond.yp = yp0 + auxPYp
        partTraj.partInitCond.gamma = gamma0*(1 + relEnSpr*randAr[4])

        partTraj = srwl.CalcPartTraj(partTraj, _mag, [1])

        if treatCT is False:
            #Search for longitudinal position:
            iz = ind0
            dZprev = 1.e+23
            newIndFound = False
            while(iz >= 0):
                if(iz < np_in_trj_mi_1):
                    if((partTraj.arZ[iz] <= _z_or_ct) and (_z_or_ct < partTraj.arZ[iz + 1])):
                        ind0 = iz
                        newIndFound = True
                        break
                dZ = abs(_z_or_ct - partTraj.arZ[iz])
                if(dZ > dZprev): break

                dZprev = dZ
                iz -= 1

            if newIndFound is False:
                iz = ind0
                dZprev = 1.e+23
                while(iz < np_in_trj_mi_1):
                    if((partTraj.arZ[iz] <= _z_or_ct) and (_z_or_ct < partTraj.arZ[iz + 1])):
                        ind0 = iz
                        newIndFound = True
                        break
                    dZ = abs(_z_or_ct - partTraj.arZ[iz])
                    if(dZ > dZprev): break

                    dZprev = dZ
                    iz += 1

            if(newIndFound is False):
                raise Exception("Trajectory point corresponding to given longitudinal position was not found")

            r0 = (_z_or_ct - partTraj.arZ[ind0])/(partTraj.arZ[ind0 + 1] - partTraj.arZ[ind0])

        x1 = partTraj.arX[ind0]; x2 = partTraj.arX[ind0 + 1]
        x = x1 + r0*(x2 - x1)
        xp1 = partTraj.arXp[ind0]; xp2 = partTraj.arXp[ind0 + 1]
        xp = xp1 + r0*(xp2 - xp1)
        y1 = partTraj.arY[ind0]; y2 = partTraj.arY[ind0 + 1]
        y = y1 + r0*(y2 - y1)
        yp1 = partTraj.arYp[ind0]; yp2 = partTraj.arYp[ind0 + 1]
        yp = yp1 + r0*(yp2 - yp1)

        xSum += x; ySum += y
        xpSum += xp; ypSum += yp
        xe2Sum += x*x; ye2Sum += y*y
        xpe2Sum += xp*xp; ype2Sum += yp*yp
        xxpSum += x*xp; yypSum += y*yp
    
    invN = 1./_npart
    xAvg = xSum*invN; yAvg = ySum*invN
    xpAvg = xpSum*invN; ypAvg = ypSum*invN
    xe2Avg = xe2Sum*invN; ye2Avg = ye2Sum*invN
    xpe2Avg = xpe2Sum*invN; ype2Avg = ype2Sum*invN
    xxpAvg = xxpSum*invN; yypAvg = yypSum*invN

    sigXe2 = xe2Avg - xAvg*xAvg; sigYe2 = ye2Avg - yAvg*yAvg
    mXXp = xxpAvg - xAvg*xpAvg; mYYp = yypAvg - yAvg*ypAvg
    sigXpe2 = xpe2Avg - xpAvg*xpAvg; sigYpe2 = ype2Avg - ypAvg*ypAvg

    return [[xAvg, xpAvg, sigXe2, mXXp, sigXpe2], [yAvg, ypAvg, sigYe2, mYYp, sigYpe2]]

#****************************************************************************
##def srwl_mag_extrap_grad_curv(_mag_curv_mid, _cen_trj_data, _xi, _xf, _nx, _yi, _yf, _ny, _zi, _zf, _nz, _grad_mult=1):
##    """
##    Extrapolates magnetic field on rectangulat 3D mesh in and off horizontal mid-plane based on field gradient in the mid-plane (e.g. for dipole with gradient or for a guad).
##    :param _mag_curv_mid: input tabulated 3D magnetic field on 2D mesh vs x and z in Natural frame in the horizontal mid-plane in a container (object of SRWLMagFldC type)
##    :param _cen_trj_data: central trajectory data defining the Natural frame of the input magnetic field (_cen_trj_data[0] is ct, _cen_trj_data[1] is x. _cen_trj_data[2] is z)
##    :param _xi: horizontal initial position of the final 3D mesh [m]
##    :param _xf: horizontal final position of the final 3D mesh [m]
##    :param _nx: number of points vs horizontal position in the final 3D mesh
##    :param _yi: vertical initial position of the final 3D mesh [m]
##    :param _yf: vertical initial position of the final 3D mesh [m]
##    :param _ny: number of points vs vertical position in the final 3D mesh
##    :param _zi: longitudinal initial position of the final 3D mesh [m]
##    :param _zf: longitudinal final position of the final 3D mesh [m]
##    :param _nz: number of points vs longitudinal position in the final 3D mesh
##    :param _grad_mult: a number the gradient has to be multiplied by
##    :returns: exrapolated 3D magnetic field structure in a container (object of SRWLMagFldC type)
##    """
##
##    if((_mag_curv_mid is None) or (_mag_curv_mid.arMagFld is None) or (len(_mag_curv_mid.arMagFld) < 1) or (_cen_trj_data is None) or (len(_cen_trj_data) < 3)
##       or (_nx < 1) or (_ny < 1) or (_nz < 1)):
##        raise Exception("Incorrect input parameters for magnetic field to be extrapolated")
##    #print('In srwl_fld_extrap_grad_curv')
##
##    fld3d_curv_mid = _mag_curv_mid.arMagFld[0]
##    #arBx_curv_mid = fld3d_curv_mid.arBx
##    arBy_curv_mid = fld3d_curv_mid.arBy
##    arBz_curv_mid = fld3d_curv_mid.arBz
