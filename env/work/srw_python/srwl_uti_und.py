#############################################################################
# SRWLib for Python: Undulator Utilities v 0.0
#############################################################################

from srwlib import *
#from copy import *

#****************************************************************************
def srwl_und_cor_fld_int(_mag3d, _dist_bw_kicks, _rms_len_kicks=0.05, _zc=0, _zcMesh=0, _zRange=0, _dupl=False):
    """
    Compensates 1st and 2nd integrals of undulator field by adding "kicks" before and after the undulator.
    :param _mag3d: 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _dist_bw_kicks: distance between kicks (in longitudinal direction) [m]
    :param _rms_len_kicks: RMS kick length [m]
    :param _zc: center position for the kick set (~center of indulator) [m]
    :param _zcMesh: longitudinal center position for the magnetic field mesh [m]
    :param _zRange: range of magnetic field to be taken into account [m]
    :param _dupl: duplicate the magnetic field object or not
    :returns: undated magnetic 3d field structure
    """

    if(_mag3d.nz <= 1): return _mag3d

    nz = _mag3d.nz
    zStep = _mag3d.rz/(nz - 1)

    #print(0.5*_dist_bw_kicks, _rms_len_kicks, _zc)

    halfDistBwKicks = 0.5*_dist_bw_kicks
    zIn = _zc - halfDistBwKicks #Position of Input Kick
    zOut = _zc + halfDistBwKicks #Position of Output Kick

    #zE = 0.5*_mag3d.rz #? End position of magnetic field
    #zB = -zE #? Start position of magnetic field
    halfFullRangeZ = 0.5*_mag3d.rz
    zB = _zcMesh - halfFullRangeZ #? End position of magnetic field
    zE = _zcMesh + halfFullRangeZ #? Start position of magnetic field

    print('zIn=', zIn, 'zOut=', zOut, 'zc=', _zc)
    print('zB=', zB, 'zE=', zE)

    zBr = zB
    zEr = zE
    if(_zRange > 0.):
        halfRange = 0.5*_zRange
        zBr = _zcMesh - halfRange
        zEr = _zcMesh + halfRange
    
    invDistBwKicks = 1./_dist_bw_kicks

    arBxOrig = _mag3d.arBx
    arByOrig = _mag3d.arBy
    manyTrPos = False
    if((_mag3d.nx > 1) or (_mag3d.ny > 1)):
        manyTrPos = True
        arBxOrig = array('d', [0]*nz) 
        arByOrig = array('d', [0]*nz) 

    halfRangeZ = 4.*_rms_len_kicks
    izInB = round((zIn - halfRangeZ - zB)/zStep)
    if(izInB < 0): izInB = 0
    izInE = round((zIn + halfRangeZ - zB)/zStep)
    if(izInE >= nz): izInE = nz - 1
    izOutB = round((zOut - halfRangeZ - zB)/zStep)
    if(izOutB < 0): izOutB = 0
    izOutE = round((zOut + halfRangeZ - zB)/zStep)
    if(izOutE >= nz): izOutE = nz - 1

    izBr = round((zBr - zB)/zStep)
    if(izBr < 0): izBr = 0
    izEr = round((zEr - zB)/zStep)
    if(izEr >= nz): izEr = nz - 1

    resMag3d = _mag3d
    if(_dupl == True): resMag3d = deepcopy(_mag3d)
    arBxRes = resMag3d.arBx
    arByRes = resMag3d.arBy

    #print(izInB, izInE, izOutB, izOutE)
    #print(izInB, izInE, izOutB, izOutE)
    print(izBr, izEr)
    
    perY = _mag3d.nx
    perZ = perY*_mag3d.ny
    for iy in range(_mag3d.ny):
        iy_perY = iy*perY
        for ix in range(_mag3d.nx):
            if(manyTrPos == True):
                arBxTot = _mag3d.arBx
                arByTot = _mag3d.arBy
                for iz in range(nz):
                    ofst = ix + iy_perY + iz*perZ
                    if(arBxTot != None):
                        arBxOrig = arBxTot[ofst]
                    if(arByTot != None):
                        arByOrig = arByTot[ofst]

            if(arBxOrig != None):
                arBxInt = uti_math.integ_array(arBxOrig, zStep, True)
                
                #auxI1X = -arBxInt[nz - 1] #[T.m]
                auxI1X_Er = arBxInt[izEr]
                auxI1X_Br = arBxInt[izBr]
                auxI1X = -(auxI1X_Er - auxI1X_Br) #[T.m]
                
                arBxInt = uti_math.integ_array(arBxInt, zStep, True)
                #auxI2X = -arBxInt[nz - 1] #[T.m^2]
                auxI2X = -(arBxInt[izEr] - arBxInt[izBr] - auxI1X_Br*(izEr - izBr)*zStep) #[T.m^2]
                
                kickI1Xin = (auxI2X - auxI1X*(zE - zOut))*invDistBwKicks
                kickI1Xout = (auxI1X - kickI1Xin)

                #print(auxI1X, 'T.m')
                #print(kickI1Xin*(zE - zIn), kickI1Xout*(zE - zOut))
                #I2test = kickI1Xin*(zE - zIn) + kickI1Xout*(zE - zOut)
                #print(auxI2X, I2test)

                B0kickX = kickI1Xin*0.3989422804/_rms_len_kicks
                z = zB + zStep*izInB
                for iz in range(izInB, izInE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zIn)/_rms_len_kicks
                    arBxRes[ofst] += B0kickX*exp(-0.5*t*t)
                    z += zStep
                B0kickX = kickI1Xout*0.3989422804/_rms_len_kicks
                z = zB + zStep*izOutB
                for iz in range(izOutB, izOutE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zOut)/_rms_len_kicks
                    arBxRes[ofst] += B0kickX*exp(-0.5*t*t)
                    z += zStep
                
            if(arByOrig != None):
                arByInt = uti_math.integ_array(arByOrig, zStep, True)
                
                #auxI1Y = -arByInt[nz - 1] #[T.m]
                auxI1Y_Er = arByInt[izEr]
                auxI1Y_Br = arByInt[izBr]
                auxI1Y = -(auxI1Y_Er - auxI1Y_Br) #[T.m]
                
                arByInt = uti_math.integ_array(arByInt, zStep, True)
                #auxI2Y = -arByInt[nz - 1] #[T.m^2]
                auxI2Y = -(arByInt[izEr] - arByInt[izBr] - auxI1Y_Br*(izEr - izBr)*zStep) #[T.m^2]

                kickI1Yin = (auxI2Y - auxI1Y*(zE - zOut))*invDistBwKicks
                kickI1Yout = (auxI1Y - kickI1Yin)

                print(auxI1Y, 'T.m')
                #I2test = kickI1Yin*(zE - zIn) + kickI1Yout*(zE - zOut)
                #print(auxI2Y, I2test)

                B0kickY = kickI1Yin*0.3989422804/_rms_len_kicks
                z = zB + zStep*izInB
                for iz in range(izInB, izInE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zIn)/_rms_len_kicks
                    arByRes[ofst] += B0kickY*exp(-0.5*t*t)
                    z += zStep
                B0kickY = kickI1Yout*0.3989422804/_rms_len_kicks
                z = zB + zStep*izOutB
                for iz in range(izOutB, izOutE + 1):
                    ofst = ix + iy_perY + iz*perZ
                    t = (z - zOut)/_rms_len_kicks
                    arByRes[ofst] += B0kickY*exp(-0.5*t*t)
                    z += zStep
    return resMag3d

#****************************************************************************
def srwl_und_fld_add_const(_mag3d, _zcMesh=0, _zc=0, _zRange=0, _bx=0, _by=0, _bz=0, _dupl=False):
    """
    Adds constant magnetic field within given longitudinal range.
    :param _mag3d: 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _zcMesh: longitudinal center position for the magnetic field mesh [m]
    :param _zc: center position for the dipole field range [m]
    :param _zRange: range of magnetic field over which to add constant field [m]
    :param _bx: horizontal field component to add [T]
    :param _by: vertical field component to add [T]
    :param _bz: longitudinal field component to add [T]
    :param _dupl: duplicate the magnetic field object or not
    :returns: undated magnetic 3d field structure
    """
    
    if(_mag3d.nz <= 1): return _mag3d

    nz = _mag3d.nz
    zStep = _mag3d.rz/(nz - 1)
    inv_zStep = 0. if(zStep == 0.) else 1./zStep

    izAddB = 0
    izAddE = nz - 1

    if(_zRange > 0):
        halfFullRangeZ = 0.5*_mag3d.rz
        zB = _zcMesh - halfFullRangeZ #? End position of magnetic field
        zE = _zcMesh + halfFullRangeZ #? Start position of magnetic field

        halfRange = 0.5*_zRange
        zAddB = _zc - halfRange
        zAddE = _zc + halfRange

        izAddTestB = int(round((zAddB - zB)*inv_zStep))
        izAddTestE = int(round((zAddE - zB)*inv_zStep))

        if(izAddTestB <= izAddTestE):
            if((izAddB < izAddTestB) and (izAddTestB <= izAddE)): izAddB = izAddTestB
            if((izAddB <= izAddTestE) and (izAddTestE < izAddE)): izAddE = izAddTestE

    resMag3d = _mag3d
    if(_dupl == True): resMag3d = deepcopy(_mag3d)

    arBX = resMag3d.arBx
    arBY = resMag3d.arBy
    arBZ = resMag3d.arBz

    bxIsSet = False
    if(arBX != None):
        if(len(arBX) > 0): bxIsSet = True 
    byIsSet = False
    if(arBY != None):
        if(len(arBY) > 0): byIsSet = True
    bzIsSet = False
    if(arBZ != None):
        if(len(arBZ) > 0): bzIsSet = True
    
    for iz in range(izAddB, izAddE + 1):
        if(bxIsSet == True): arBX[iz] += _bx
        if(byIsSet == True): arBY[iz] += _by
        if(bzIsSet == True): arBZ[iz] += _bz

    return resMag3d

#****************************************************************************
def srwl_und_cut_fld(_mag3d, _res_rz, _zc=None, _dupl=False):
    """
    Cuts (truncates) undulator magnetic field. Assumes equidistant mesh.
    :param _mag3d: 3D magnetic field of undulator (object of SRWLMagFld3D type)
    :param _res_rz: range vs longitudinal position to produce [m]
    :param _zc: center position for the kick set (~center of indulator) [m]
    :param _dupl: duplicate the magnetic field object or not
    :returns: undated magnetic 3d field structure
    """

    if(_mag3d.nz <= 1): return _mag3d

    nx = _mag3d.nx
    ny = _mag3d.ny
    nz = _mag3d.nz
    zStep = _mag3d.rz/(nz - 1)
    zRange = zStep*(nz - 1)

    if((_res_rz <= 0) or (_res_rz >= zRange)):
        if(_dupl == False): return _mag3d
        else: return deepcopy(_mag3d)

    zBegOrig = -0.5*zRange
    zEndOrig = -zBegOrig

    perZ = nx*ny

    arBxOrig = _mag3d.arBx
    arByOrig = _mag3d.arBy
    arBzOrig = _mag3d.arBz

    if(_zc == None):
        #Find Max. field value and position
        Be2max = 0.
        i = 0
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    Bx = 0 if arBxOrig == None else arBxOrig[i]
                    By = 0 if arByOrig == None else arByOrig[i]
                    Bz = 0 if arBzOrig == None else arBzOrig[i]
                    curBe2 = Bx*Bx + By*By + Bz*Bz
                    if(Be2max < curBe2): Be2max = curBe2
                    i += 1

        #Find Longitudinal Center position of the field 
        Bthresh = sqrt(Be2max)*0.1 #to steer
        #print('Bthresh=', Bthresh)

        i = 0; izStart = 0
        wasBreak = False
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    Bx = 0 if arBxOrig == None else arBxOrig[i]
                    By = 0 if arByOrig == None else arByOrig[i]
                    Bz = 0 if arBzOrig == None else arBzOrig[i]
                    curB = sqrt(Bx*Bx + By*By + Bz*Bz)
                    #print('i=', i, 'curB=', curB)
                    
                    if(curB >= Bthresh):
                        izStart = iz
                        wasBreak = True
                        break
                    i += 1
                if(wasBreak): break
            if(wasBreak): break

        izEnd = nz - 1
        wasBreak = False
        for iz in range(nz - 1, -1, -1):
            iz_perZ = iz*perZ
            for iy in range(ny):
                iy_nx_p_iz_perZ = iy*nx + iz_perZ
                for ix in range(nx):
                    i = ix + iy_nx_p_iz_perZ
                    Bx = 0 if arBxOrig == None else arBxOrig[i]
                    By = 0 if arByOrig == None else arByOrig[i]
                    Bz = 0 if arBzOrig == None else arBzOrig[i]
                    curB = sqrt(Bx*Bx + By*By + Bz*Bz)
                    if(curB >= Bthresh):
                        izEnd = iz
                        wasBreak = True
                        break
                if(wasBreak): break
            if(wasBreak): break

        print('izStart=', izStart, 'izEnd=', izEnd)

        zStart = zBegOrig + izStart*zStep
        zEnd = zBegOrig + izEnd*zStep
        _zc = 0.5*(zStart + zEnd)

    halfResRangeZ = 0.5*_res_rz
    zStartRes = _zc - halfResRangeZ
    if(zStartRes < zBegOrig): zStartRes = zBegOrig
    zEndRes = _zc + halfResRangeZ
    if(zEndRes > zEndOrig): zEndRes = zEndOrig

    zRangeRes = zEndRes - zStartRes

    izStartRes = int(round((zStartRes - zBegOrig)/zStep))
    if(izStartRes < 0): izStartRes = 0

    izEndRes = int(round((zEndRes - zBegOrig)/zStep))
    if(izEndRes >= nz): izEndRes = nz - 1

    if(izEndRes < izStartRes): izEndRes = izStartRes

    nzRes = izEndRes - izStartRes + 1
    nTot = perZ*nzRes
    auxList = [0]*nTot
    arBxRes = None if(arBxOrig == None) else array('d', auxList)
    arByRes = None if(arByOrig == None) else array('d', auxList)
    arBzRes = None if(arBzOrig == None) else array('d', auxList)

    for iz in range(izStartRes, izEndRes + 1):
        iz_perZ = iz*perZ
        iz0_perZ = (iz - izStartRes)*perZ
        for iy in range(ny):
            iy_nx_p_iz_perZ = iy*nx + iz_perZ
            iy_nx_p_iz0_perZ = iy*nx + iz0_perZ
            
            for ix in range(nx):
                i = ix + iy_nx_p_iz_perZ
                i0 = ix + iy_nx_p_iz0_perZ
                if(arBxOrig != None): arBxRes[i0] = arBxOrig[i]
                if(arByOrig != None): arByRes[i0] = arByOrig[i]
                if(arBzOrig != None): arBzRes[i0] = arBzOrig[i]
                
    if(_dupl==True):
        return SRWLMagFld3D(_arBx=arBxRes, _arBy=arByRes, _arBz=arBzRes, _nx=nx, _ny=ny, _nz=nzRes,
                            _rx=_mag3d.rx, _ry=_mag3d.ry, _rz=zRangeRes, _nRep=_mag3d.nRep, _interp=_mag3d.interp)
    else:
        _mag3d.arBx = arBxRes; _mag3d.arBy = arByRes; _mag3d.arBz = arBzRes; _mag3d.nz = nzRes; _mag3d.rz = zRangeRes

    return _mag3d

#****************************************************************************
##def srwl_und_max_flux_spec(_e_beam, _mesh, _per, _len, _bx_vs_par, _by_vs_par, [_par1_start, _par1_step], [_par1_start, _par1_step], _par1_to_use, _par2_step_to_use, _min_pol_rate):
##    """
##    Calculate Stokes Parameters of Emitted (and Propagated, if beamline is defined) Partially-Coherent SR
##    :param _e_beam: Finite-Emittance e-beam (SRWLPartBeam type)
##    :param _mesh: mesh vs photon energy, horizontal and vertical positions (SRWLRadMesh type) on which initial SR should be calculated
##    :param _per: undulator period
##    :param _len: undulator length
##    :param _n_part_tot: total number of "macro-electrons" to be used in the calculation
##    :param _n_part_avg_proc: number of "macro-electrons" to be used in calculation at each "slave" before sending Stokes data to "master" (effective if the calculation is run via MPI)
##    :param _n_save_per: periodicity of saving intermediate average Stokes data to file by master process
##    :param _file_path: path to file for saving intermediate average Stokes data by master process
##    :param _sr_samp_fact: oversampling factor for calculating of initial wavefront for subsequent propagation (effective if >0)
##    :param _opt_bl: optical beamline (container) to propagate the radiation through (SRWLOptC type)
##    :param _pres_ang: switch specifying presentation of the resulting Stokes parameters: coordinate (0) or angular (1)
##    :param _char: radiation characteristic to calculate: 0- Intensity (s0); 1- Four Stokes components; 2- Mutual Intensity Cut vs X; 3- Mutual Intensity Cut vs Y; 4- Mutual Intensity Cut vs X & Y
##    :param _x0: horizontal center position for mutual intensity calculation
##    :param _y0: vertical center position for mutual intensity calculation
##    """

#*********************************Entry point (for command-prompt calls)
if __name__=="__main__":
    import optparse
    p = optparse.OptionParser()
    p.add_option('--cor', dest='do_cor_fld_int', action='store_true', default=False, help="correct 1st and 2nd field integrals of tabulated undulator magnetic field")
    p.add_option('--conv', dest='do_conv_fld', action='store_true', default=False, help="convert magnetic measurement data file to SRW ASCII format")
    p.add_option('--ifn', dest='ifn', type='string', default='', help="input field file name")
    p.add_option('--ofn', dest='ofn', type='string', default='', help="output field file name")
    p.add_option('--otfn', dest='otfn', type='string', default='', help="output trajectory file name")
    p.add_option('--dbwk', dest='dbwk', type="float", default=2.8, metavar="NUMBER", help="distance between correcting kicks [m]")
    p.add_option('--rmsk', dest='rmsk', type="float", default=0.05, metavar="NUMBER", help="RMS length of the correcting kicks [m]")
    p.add_option('--zk', dest='zk', type="float", default=0., metavar="NUMBER", help="longitudinal position of the correcting kicks [m]")
    p.add_option('--z0t', dest='z0t', type="float", default=1.e+23, metavar="NUMBER", help="initial longitudinal position for trajectory calculation [m]")
    p.add_option('--x0t', dest='x0t', type="float", default=0., metavar="NUMBER", help="initial horizontal position for trajectory calculation [m]")
    p.add_option('--y0t', dest='y0t', type="float", default=0., metavar="NUMBER", help="initial vertical position for trajectory calculation [m]")
    p.add_option('--xp0t', dest='xp0t', type="float", default=0., metavar="NUMBER", help="initial horizontal angle for trajectory calculation [m]")
    p.add_option('--yp0t', dest='yp0t', type="float", default=0., metavar="NUMBER", help="initial vertical angle for trajectory calculation [m]")
    p.add_option('--elen', dest='elen', type="float", default=3., metavar="NUMBER", help="electron energy [GeV]")
    p.add_option('--dbx', dest='dbx', type="float", default=0., metavar="NUMBER", help="horizontal magnetic field component to add [T]")
    p.add_option('--dby', dest='dby', type="float", default=0., metavar="NUMBER", help="vertical magnetic field component to add [T]")
    p.add_option('--dbz', dest='dbz', type="float", default=0., metavar="NUMBER", help="longitudinal magnetic field component to add [T]")
    p.add_option('--dzc', dest='dzc', type="float", default=0., metavar="NUMBER", help="center longitudinal position for adding constant magnetic field [m]")
    p.add_option('--dzr', dest='dzr', type="float", default=0., metavar="NUMBER", help="range of longitudinal position within which to add constant magnetic field [m]")
    
    opt, args = p.parse_args()

    #print(opt.ifn)
    fldCnt = None
    
    if(opt.do_conv_fld): #Convertion from magnetic measurements data format
        fldDataCols = srwl_uti_read_data_cols(opt.ifn, '\t')
        arZ = fldDataCols[0]
        zNp = len(arZ)
        zSt = arZ[0]*0.001
        zFi = arZ[zNp - 1]*0.001
        fldCnt = SRWLMagFldC(SRWLMagFld3D(array('d', fldDataCols[1]), array('d', fldDataCols[2]), array('d', fldDataCols[3]), 1, 1, zNp, 0., 0., zFi - zSt), 0., 0., 0.5*(zSt + zFi))

    if(fldCnt == None): fldCnt = srwl_uti_read_mag_fld_3d(opt.ifn)
    fld3d = fldCnt.arMagFld[0]
    zc = fldCnt.arZc[0]

    if((opt.dbx != 0.) or (opt.dby != 0.) or (opt.dbz != 0.)):
        srwl_und_fld_add_const(fld3d, fldCnt.arZc[0], opt.dzc, opt.dzr, opt.dbx, opt.dby, opt.dbz)

    if(opt.do_cor_fld_int):
        zc = fldCnt.arZc[0]
        srwl_und_cor_fld_int(fld3d, opt.dbwk, opt.rmsk, opt.zk, zc)

    if(len(opt.ofn) > 0):
        fld3d.save_ascii(opt.ofn, 0., 0., zc)

    if(len(opt.otfn) > 0):
        z0 = opt.z0t if(opt.z0t != 1.e+23) else fldCnt.arZc[0] - 0.5*fld3d.rz
        elec = SRWLParticle(_x=opt.x0t, _y=opt.y0t, _z=z0, _xp=opt.xp0t, _yp=opt.yp0t, _gamma=opt.elen/0.51099890221e-03)
        traj = SRWLPrtTrj(_np=fld3d.nz, _ctStart=0, _ctEnd=fld3d.rz, _partInitCond=elec)
        #print(fld3d.arBy)
        srwl.CalcPartTraj(traj, fldCnt, [1])
        traj.save_ascii(opt.otfn)
