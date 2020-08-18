#############################################################################
#This is a subset of Gianluca's and Ilya's module to read-in Genesis output files
#and to generate SRW wavefront structure
#############################################################################

#from mpi4py import MPI
#comm  = MPI.COMM_WORLD   #Creates one common object "communicator" knwoning about processes
#rank  = comm.Get_rank()  #Process number
#nproc = comm.Get_size()  #Total number of processes

from __future__ import print_function #Python 2.7 compatibility
from srwlib import *

#from xframework.adaptors.genesis import *
##################from ocelot.adaptors.genesis import *

#from prop_lines import *
import copy
import struct

#from scipy import interpolate

##################################################################################################################
############################################ ROUTINES FROM OCELOT ################################################
####                                                                                                          ####
#### These are to read the Genesis Output file; necessary here in order to make the module 'standalone' ##########
#######                                          i.e. without Ocelot                                    ##########
##################################################################################################################              

#import numpy as np
#from copy import copy, deepcopy

class GenesisOutput:
    
    def __init__(self):
        self.z = []
        self.I = []
        self.n = []
        self.zSlice = []
        self.E = []
        self.aw = []
        self.qfld = []
        
        self.sliceKeys = []
        self.sliceValues = {}
        
        self.parameters = {}
    
    def __call__(self, name):
        '''
        if name not in self.__dict__.keys():
            return 0
        else:
            return self.__dict__[name]
        '''
        if name not in self.parameters.keys():
            return 0.0
        else:
            p, = self.parameters[name]
            return float(p.replace('D','E'))
        

#******************************************************************************
def uti_io_read_Genesis_output(fileName, _short=True):
#def readGenesisOutput(fileName, ):

    out = GenesisOutput()
    
    chunk = 'header'
    
    #nSlice = 0
    
    f=open(fileName,'r')
    f.readline()
    for line in f:

        #print(line)
        
        #tokens = line.replace('=', ' ').split()
        tokens = line.strip().split()

        if len(tokens) < 1:
            #chunk = 'none'
            continue
        
        if tokens == ['z[m]', 'aw', 'qfld']:

            if(_short): return out #OC
            
            chunk = 'optics'
            #print 'reading optics '
            print('reading optics')
            continue
        
        if tokens[0] == 'Input':
            chunk = 'input'
            #print 'reading input parameters'
            #print('reading input parameters')
            continue
        
        #********** output: slice    10
        if tokens[0] == '**********':
            #print 'slice:', tokens[3]
            chunk = 'slices'
            nSlice = int(tokens[3])
         
        if tokens[0] == 'power':
            chunk = 'slice'
            out.sliceKeys = copy.copy(tokens)
            #out.sliceValues[nSlice] = copy.copy(tokens)
            out.sliceValues[nSlice] = {}
            for i in range(0,len(tokens)):
                out.sliceValues[nSlice][out.sliceKeys[i]] = []
            #print 'reading slices'
            #print out.sliceKeys
            #print out.sliceValues[nSlice]
            continue
            
        if chunk == 'optics':
            z,aw,qfld = map(float,tokens)
            out.z.append(z)
            out.aw.append(aw)
            out.qfld.append(qfld)
            
        if chunk == 'input':
            #tokens=line.replace('=','').strip().split()
            tokens=line.replace('=',' ').strip().split() #OC: fix to ensure correct parsing of tokens like 'nslice=10932'
            out.parameters[tokens[0]] = tokens[1:]
            #print('input:', tokens)

        if chunk == 'slice':

            vals = map(float,tokens)
            if sys.hexversion >= 0x030000F0: vals = list(vals) #OC: for compatibility with Py3
            
            #print vals
            for i in range(0,len(vals)):
                out.sliceValues[nSlice][out.sliceKeys[i]].append(vals[i])
            
            #out.zSlice.append(vals[2])
            #out.aw.append(aw)
            #out.qfld.append(qfld)

        if chunk == 'slices':
            if len(tokens) == 2 and tokens[1]=='current':
                #print tokens[1]
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)

    #out.nSlices = int(out('nslice'))
    #out.nSlices = len(out.sliceValues) 
    #out.nZ = len(out.sliceValues[1][out.sliceKeys[0]])

    #print 'nSlice', out.nSlices
    #print 'nZ', out.nZ
    #print('nSlice', out.nSlices)
    #print('nZ', out.nZ)

    return out

#******************************************************************************
#def uti_io_read_Genesis_rad(fileName='simulation.gout.dfl', npoints=51, slice_start=0, slice_end = -1, idx=None):
def uti_io_read_Genesis_rad(fileName='simulation.gout.dfl', _ncar=51, _mult=1., _slice_skip_per=0, _slice_start=0, _slice_end=-1, _align_out='txy'):
    
    #def read_in_chunks(_file, _size=1024):
    #    while True:
    #        data = _file.read(_size)
    #        if not data:
    #            break
    #        yield data

    #def read_in_data_block(_file, _ofst=0, _size=1024):
    #    if(_ofst > 0): _file.seek(_ofst)
    #    return _file.read(_size)

    ncarE2 = _ncar*_ncar
    slice_size = int(ncarE2*16)

    f = open(fileName,'rb')
    
    f.seek(0,2)
    #total_data_len = f.tell() / 8 / 2
    total_data_len = f.tell() / 16 # = _ncar*_ncar*n_slice

    nSliceTot = int(round(total_data_len/ncarE2))
    iSliceEnd = int(_slice_end)
    if(_slice_end < 0): iSliceEnd = nSliceTot - 1

    nSliceToRead = int((iSliceEnd + 1 - _slice_start)/(_slice_skip_per + 1))
    nResVal = nSliceToRead*ncarE2*2
    
    #print('starting array allocation...')
    resFldData = srwl_uti_array_alloc('f', nResVal)
    #print('done')

    #if _slice_start > 0:
    #    f.seek(slice_size*slice_start,0)

    #if slice_end > 0:
    #    data_size = (slice_end - slice_start) * slice_size
    #else:
    #    data_size = total_data_len * 16 - (slice_start) * slice_size

    #print('slice size = ', slice_size)
    #print('data size = ', data_size)
    #print('total_data_len = ', total_data_len)
    
    #ncar = int(np.sqrt((slice_size/16)))
    #print('ncar=', ncar)
    
    #n_slices = total_data_len / (slice_size/16)
    
    #if idx == None:
        #slices = np.zeros([n_slices,ncar,ncar], dtype=complex)
        #slices = array('f', [0]*(n_slices*ncar*ncar*2))
        #slices = srwl_uti_array_alloc('f', (total_data_len*2))
    #else:
        #slices = np.zeros([len(idx),ncar,ncar], dtype=complex)
        #slices = array('f', [0]*(len(idx)*ncar*ncar*2))
        #slices = srwl_uti_array_alloc('f', (len(idx)*ncarE2*2))

    perT = 2
    perX = nSliceToRead*perT
    perY = perX*_ncar

    if((_align_out == 'xyt') or (_align_out == 'xye')):
        perX = 2
        perY = perX*_ncar
        perT = perY*_ncar
    elif((_align_out == 'yxt') or (_align_out == 'yxe')):
        perY = 2
        perX = perY*_ncar
        perT = perX*_ncar
    #else...

    nSliceToRead_mi1 = nSliceToRead - 1

    ofstFile = _slice_start*slice_size
    ofstFileStep = _slice_skip_per*slice_size
    f.seek(ofstFile, 0)

    for iSlice in range(nSliceToRead):
        piece = f.read(slice_size)
        iSlice_perT = iSlice*perT
        i = 0
        for iy in range(_ncar):
            iy_perY = iy*perY
            for ix in range(_ncar):
                ofstRes = int(iSlice_perT + ix*perX + iy_perY)
                i16 = i*16
                i16p8 = i16 + 8
                resFldData[ofstRes] = _mult*(struct.unpack('d', piece[i16 : i16p8])[0])
                resFldData[ofstRes + 1] = _mult*(struct.unpack('d', piece[i16p8 : i16p8 + 8])[0])
                i += 1
        if(_slice_skip_per > 0):
            if(iSlice < nSliceToRead_mi1): f.seek(ofstFileStep, 1) #move to next read position

        #print('slices read:', iSlice)
    return resFldData

#******************************************************************************
def uti_io_set_wfr_from_Genesis_files(_fpGenOut, _fpGenRad, _pol='h', _slice_skip_per=0, _slice_start=0, _slice_end=-1):
    """Sets up SRW wavefront from Genesis output
    :param _fpGenOut: file path to Genesis Output file
    :param _fpGenRad: file path to Genesis Radiation file
    :param _pol: polarization ('h' for horizontal, 'v' for vertical,...)
    :param _slice_skip_per: number of consequtive slice to skip after reading-in one slice
    :param _slice_start: first slice number (in the Genesis file) to read in
    :param _slice_end: last slice number (in the Genesis file) to read in
    """

    outGenesis = uti_io_read_Genesis_output(_fpGenOut)
    
    ncar = int(outGenesis('ncar'))
    zrayl = outGenesis('zrayl')
    rmax0 = outGenesis('rmax0')
    xlamds = outGenesis('xlamds')
    rxbeam = outGenesis('rxbeam')
    rybeam = outGenesis('rybeam')
    zsep = outGenesis('zsep')
    ntail = outGenesis('ntail')
    itdp = int(outGenesis('itdp'))
    xlamd = outGenesis('xlamd')

    pi = 3.14159265358979
    lightSpeed = 2.9979245812E+08 #Speed of Light [m/c]
    light_eV_mu = 1.23984186
    #vacimp = 376.7303135 #Vacuum impedance in Ohms
    #eev = 510998.902 #Energy units (mc^2) in eV

    photEn_xlamds = light_eV_mu*1.e-06/xlamds

    w0 = sqrt(zrayl*xlamds/pi)
    sigrB = sqrt(rxbeam*rxbeam + rybeam*rybeam)
    #sigrF = w0
    dGrid = 0.5*rmax0*(sigrB + w0) #Mesh from -grid to +grid; This is the genesis input parameter dgrid/2

    #This defines a constant to be used for conversion (multiplication) of GENESIS Electric field
    #to SRW TD Electric Field in units of sqrt(W/mm^2)
    #xkper0 = 2*pi/xlamd #undulator "wavenumber"
    #xks = 2*pi/xlamds
    #convE_Gen2SRW = (1e-03)*xkper0*xkper0*eev/(xks*sqrt(vacimp))

    xStep = 2.*dGrid/(ncar - 1)
    convE_Gen2SRW = sqrt(1.e-06/(xStep*xStep)) #? See Genesis 1.3 manual, page 41
    #convE_Gen2SRW = 1. #?
    #print('convE_Gen2SRW=', convE_Gen2SRW)

    arEx = uti_io_read_Genesis_rad(_fpGenRad, ncar, _mult=convE_Gen2SRW, _slice_skip_per=_slice_skip_per)
    lenArE = len(arEx)
    nt = int(round(lenArE/(ncar*ncar*2)))
    arEy = None
    if(_pol == 'v'): #to treat circular polarization!?
        arEy = arEx
        arEx = None
    if(arEx == None): arEx = srwl_uti_array_alloc('f', lenArE)
    if(arEy == None): arEy = srwl_uti_array_alloc('f', lenArE)

    tStepOrig = xlamds*zsep/lightSpeed;
    tStep = tStepOrig*(_slice_skip_per + 1)
    tStartOrig = xlamds*zsep*ntail/lightSpeed;
    tStart = tStartOrig + _slice_start*tStepOrig
    tEnd = tStart + tStep*(nt - 1)

    presFT = 1 #presentation/domain: 0- frequency (photon energy), 1- time
    
    if(itdp == 0):
        tStart = photEn_xlamds #to treat harmonics!!!
        tStep = 0
        presFT = 0

    wfr = SRWLWfr(arEx, arEy, _typeE='f',
                  _eStart=tStart, _eFin=tEnd, _ne=nt,
                  _xStart=-dGrid, _xFin=dGrid, _nx=ncar,
                  _yStart=-dGrid, _yFin=dGrid, _ny=ncar, _zStart=0)
    wfr.avgPhotEn = photEn_xlamds
    wfr.presFT = presFT
    wfr.unitElFld = 2 #electric field units are sqrt(J/eV/mm^2) or sqrt(W/mm^2), depending on representation (freq. or time)

    return wfr


'''    
    n = 0
    #id = 0
    iss = 0
    
    numSlicesToReadAtOnce = 10 #to tune

    #for piece in read_in_chunks(f, size  = slice_size):
    for piece in read_in_chunks(f, size = numSlicesToReadAtOnce*slice_size):

        print('slices read:', iss)

        for iSliceRead in range(numSlicesToReadAtOnce):

            if (idx == None) or (n in idx):
                #print 'reading', n
                #if ( len(piece) / 16 != ncar**2):
                    #print 'warning, wrong slice size'
                    #print('warning, wrong slice size')
        
                #for i in xrange(len(piece) / 16 ):
                for i in range(ncarE2): #OC: for compatibility with Py3
                
                    i2 = i % _ncar
                    i1 = int(i / _ncar)
                    #print struct.unpack('d',piece[16*i:16*i+8])
                    #slices[id,i1,i2] = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0]
                    ofst = iss + i1*perX + i2*perY
                    slices[ofst] = struct.unpack('d', piece[16*i:16*i+8])[0]
                    slices[ofst + 1] = struct.unpack('d', piece[16*i+8:16*(i+1)])[0]

                    #if(slices[ofst] != 0.): print(slices[ofst], slices[ofst + 1])

                #print('slice #', iss, ' read')

                #id += 1
                iss += 1
            n += 1
    
    #print('read slices: ', slices.shape)
    return slices


def readRadiationFile_mpi(comm=None, fileName='simulation.gout.dfl', npoints=51):
    """
    not advisable to be used with very small n_proc due to memory overhead ~ file_size / n_proc  
    """
    from mpi4py import MPI
    
    def read_in_chunks(file, size=1024):
        while True:
            data = file.read(size)
            if not data:
                break
            yield data

    
    rank = comm.Get_rank()
    nproc = comm.Get_size()
        
    f = open(fileName,'rb')
    f.seek(0,2)
    total_data_len = f.tell() / 8 / 2
    f.seek(0,0)
    slice_size = int(2.0*npoints*npoints * 8.0)
    n_slices = int(total_data_len / (npoints**2))

    ncar = int(np.sqrt((slice_size/16)))
    
    local_data_len = int(n_slices / nproc)
    
    n_extra = n_slices - local_data_len * nproc
    
    tmp_buf = np.zeros([local_data_len,ncar,ncar], dtype=complex)

    
    if rank == 0:
        slice_start  = rank * local_data_len
        slices = np.zeros([n_slices,ncar,ncar], dtype=complex)
        slices_to_read = local_data_len + n_extra
    else:
        slice_start = rank * local_data_len + n_extra
        slices = []
        slices_to_read = local_data_len
        
    n = 0
    
    f.seek(slice_start*slice_size, 0)

    print 'rank', rank, ' reading', slice_start, slices_to_read, n_extra

    for piece in read_in_chunks(f, size  = slice_size):
                
        if n >= slices_to_read :
            break
        
        if ( len(piece) / 16 != ncar**2):
            print 'warning, wrong slice size'
    
        for i in xrange(len(piece) / 16 ):
            i2 = i % ncar
            i1 = int(i / ncar)
            if rank == 0:
                #print n, n_extra
                v = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0]
                slices[n,i1,i2] = v
                if n - n_extra >= 0:
                    tmp_buf[n-n_extra,i1,i2] = v
            else:
                tmp_buf[n,i1,i2] = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0] 
        n += 1
    #print rank, 'tmp_buf=', tmp_buf
    comm.Gather([tmp_buf,  MPI.COMPLEX], [slices[n_extra:], MPI.COMPLEX])
    
    return slices


#################################################################################################
####################################   End routines from Ocelot #################################
#################################################################################################

#################################################################################################
####################             Routines NOT used in prop_par.py           #####################
####################          but useful for near-future extensions         #####################
#################################################################################################

def fluence(field,Xdim,Ydim,nx,ny,Ntotph):

    Dx = Xdim/nx
    Dy = Ydim/ny
    print('Dx = ',Dx,' m')
    print('Dy = ',Dy,' m')
    integral = 0    
    maxintens = field[0]*field[0] + field[1]*field[1]    
    for j in np.arange(0,len(field),2):        
        intens = field[j]*field[j] + field[j+1]*field[j+1]
        integral = integral + intens * Dx * Dy
        if maxintens < intens: maxintens = intens        
    result = Ntotph*maxintens/integral
    print('max intensity = ',maxintens,' A.U.')
    print('integrated intensity = ',integral,' A.U.')
    print('Xdim = ',Xdim,' m')
    print('Ydim = ',Ydim,' m')
                        
    return result


def fluenceI(Inte,Xdim,Ydim,nx,ny,Ntotph):

    Dx = Xdim/nx
    Dy = Ydim/ny
    print('Dx = ',Dx,' m')
    print('Dy = ',Dy,' m')
    integral = 0    
    maxintens = Inte[0]    
    for j in np.arange(0,len(Inte),1):        
        intens = Inte[j]
        integral = integral + intens * Dx * Dy
        if maxintens < intens: maxintens = intens        
    result = Ntotph*maxintens/integral
    print('max intensity = ',maxintens,' A.U.')
    print('integrated intensity = ',integral,' A.U.')
    print('Xdim = ',Xdim,' m')
    print('Ydim = ',Ydim,' m')
                        
    return result



def writeres(namef, dataX, dataY):

    f1 = open(namef, 'w')
    for i in range(len(dataX)):
        f1.write('%s ' %(dataX[i]) + '%s' %(dataY[i]) +'\n')        
    f1.close()
    




def findsourceave(field, wfr0, zini, zfin, zstep, slicein, slicefin):
    
    reg_arzi = []
    wfrini = copy.deepcopy(wfr0)
    NX = wfrini.mesh.nx
    NY = wfrini.mesh.ny
    
    #print('nx0',nx0)
    #print('ny0',ny0)
    R_Xs0=[]
    R_Xf0=[]
    R_nx0=[]
    R_Ys0=[]
    R_Yf0=[]
    R_ny0=[]
    
    for qsl in range(slicein,slicefin):

        print('Finding the source. Processing slice nr ',qsl-slicein+1,' out of ',slicefin-slicein) 
        sl = field[qsl]
        count = 0
        for j in range(NY):
            for k in range(NX):
                wfrini.arEx[count] = np.float(np.real(sl[j][k]))
                count = count + 1
                wfrini.arEx[count] = np.float(np.imag(sl[j][k]))
                count = count + 1
        OUT = findsource(wfrini, zini, zfin, zstep)
        registerARZI = OUT[5]
        registerMESH = OUT[8]

        
        
        inte_registerARZI = []
        
        
        if qsl == slicein:
            for num in range(len(registerARZI)):
                R_Xs0.append(registerMESH[num].xStart)
                R_Xf0.append(registerMESH[num].xFin)
                R_nx0.append(registerMESH[num].nx)
                R_Ys0.append(registerMESH[num].yStart)
                R_Yf0.append(registerMESH[num].yFin)
                R_ny0.append(registerMESH[num].ny)

        for num in range(len(registerARZI)):

            Xs0 = R_Xs0[num]
            Xf0 = R_Xf0[num]
            nx0 = R_nx0[num]
            Ys0 = R_Ys0[num]
            Yf0 = R_Yf0[num]
            ny0 = R_ny0[num]
            
                  
            XMESH = np.linspace(Xs0, Xf0, nx0)
            YMESH = np.linspace(Ys0, Yf0, ny0)
            
            #print(Xs0, ' ', Xf0, ' ',nx0, ' ',ny0)                
            Xs=registerMESH[num].xStart
            Xf=registerMESH[num].xFin
            nx=registerMESH[num].nx
            Ys=registerMESH[num].yStart
            Yf=registerMESH[num].yFin
            ny=registerMESH[num].ny
            Xact = np.linspace(Xs, Xf, nx)
            Yact = np.linspace(Ys, Yf, ny)
            dataI = np.array( registerARZI[num] )
            zI    = dataI.reshape( ny, nx )
            Interdata = interpolate.interp2d(Xact, Yact, zI)
            zI_inte    = Interdata(XMESH, YMESH)
            INT_arI   = zI_inte.reshape(ny0,nx0)
            inte_registerARZI.append(np.array(INT_arI))
            #print('slice ',qsl,' len = ',len(inte_registerARZI),' with ',len(np.array(INT_arI)))


            
            #plotinte(str(qsl)+'_'+str(num)+'old_pos' ,Xs,  Xf,  nx,  Ys,  Yf,  ny,  zI)
            #plotinte(str(qsl)+'_'+str(num)+'new_pos' ,Xs0, Xf0, nx0, Ys0, Yf0, ny0, INT_arI)
            
        
        if qsl == slicein:
            reg_arzi=np.array(inte_registerARZI)
        else:
            reg_arzi = reg_arzi + np.array(inte_registerARZI)
            

    InteWaist = 0  
    registerI = []
    registerZ = OUT[1]    
    register_fwhm_X = []
    register_fwhm_Y = []


    #print(len(reg_arzi))
    #print(len(reg_arzi[0]))
    #print(len(reg_arzi[0][0]))

    #print(len(reg_arzi),'... e` 19?')
    
    InteWaist = 0
    for k in range(len(reg_arzi)):        
        intens = np.max(reg_arzi[k])
        registerI.append(intens)        
        if intens > InteWaist:
            InteWaist = intens                       
            WaistPos = OUT[1][k]
            xwaist = 1e6 * np.linspace(R_Xs0[k], R_Xf0[k], R_nx0[k] )
            ywaist = 1e6 * np.linspace(R_Ys0[k], R_Yf0[k], R_ny0[k] )
            
        x = 1e6 * np.linspace(R_Xs0[k], R_Xf0[k], R_nx0[k] )
        y = 1e6 * np.linspace(R_Ys0[k], R_Yf0[k], R_ny0[k] )
        arIpy=np.sum(reg_arzi[k],1)
        arIpyf, fwhmy=fwhm_gauss_fit(y,arIpy)        
        arIpx=np.sum(reg_arzi[k],0)
        arIpxf, fwhmx=fwhm_gauss_fit(x,arIpx)        
        register_fwhm_X.append(fwhmx)
        register_fwhm_Y.append(fwhmy) 
          
    return [registerZ,registerI,WaistPos,InteWaist,register_fwhm_X,register_fwhm_Y,reg_arzi, xwaist, ywaist]
        





def findsource(wfrini1, zini, zfin, zstep):

    Zpos = zini-zstep
    flagf    = 0
    j = 1
    registerI = []
    registerZ = []
    registerARZI = []
    register_fwhm_X = []
    register_fwhm_Y = []
    register_mesh   = []
    print('Initial wavefront position = ',wfrini1.mesh.zStart)

    while Zpos>zfin:
        
        Drift_1 = SRWLOptD(-j*zstep+zini) #Drift    
        #Wavefront Propagation Parameters with OM:
        #                    [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11]         
        ppDrift_1 =          [ 1,  1, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]   
        ppFinal =            [ 0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

        #ppDrift_1 =          [ 0,  0, 1.0,  0,  0, 1.0, 2.0, 1.0, 2.0,  0,  0,   0]   
        #ppFinal =            [ 0,  0, 1.0,  0,  0, 0.5, 0.5, 0.5, 0.5,  0,  0,   0]
        
        optBL = SRWLOptC([Drift_1],  [ppDrift_1, ppFinal])
        wfr1c = copy.deepcopy(wfrini1)                
        srwl.PropagElecField(wfr1c, optBL)
        arI = array('f', [0]*wfr1c.mesh.nx*wfr1c.mesh.ny) #"flat" array to take 2D intensity data
        srwl.CalcIntFromElecField(arI, wfr1c, 6, 0, 3, wfr1c.mesh.eStart, 0, 0) #extracts intensity        
        arP = array('d', [0]*wfr1c.mesh.nx*wfr1c.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
        srwl.CalcIntFromElecField(arP, wfr1c, 0, 4, 3, wfr1c.mesh.eStart, 0, 0) #extracts radiation phase
        register_mesh.append(wfr1c.mesh)
        dataI = np.array( arI )
        zI=dataI.reshape( wfr1c.mesh.ny, wfr1c.mesh.nx )
        registerARZI.append(zI)
        x = 1e6 * np.linspace(wfr1c.mesh.xStart, wfr1c.mesh.xFin, wfr1c.mesh.nx )
        y = 1e6 * np.linspace(wfr1c.mesh.yStart, wfr1c.mesh.yFin, wfr1c.mesh.ny )
        arIpy=np.sum(zI,1)
        arIpyf, fwhmy=fwhm_gauss_fit(y,arIpy)        
        arIpx=np.sum(zI,0)
        arIpxf, fwhmx=fwhm_gauss_fit(x,arIpx)        
        register_fwhm_X.append(fwhmx)
        register_fwhm_Y.append(fwhmy) 
        field = wfr1c.arEx
        maxintens = field[0]*field[0] + field[1]*field[1]    
        for k in np.arange(0,len(field),2):        
            intens = field[k]*field[k] + field[k+1]*field[k+1]
            if maxintens < intens: maxintens = intens
        registerI.append(maxintens)
        registerZ.append(Zpos)           
        Zpos = Zpos - zstep
        j = j+1
        #plotinte('singlePOS'+str(j) ,wfr1c.mesh.xStart, wfr1c.mesh.xFin, wfr1c.mesh.nx,  wfr1c.mesh.yStart, wfr1c.mesh.yFin, wfr1c.mesh.ny,   zI)

    InteWaist = np.max(registerI)    
    maxind = np.argmax(registerI)
    WaistPos = registerZ[maxind]

    [wfr1d, arI, arP] = Source_Prop(Len_drft = WaistPos-zini, Wfrin = wfrini1) 

    
    return [wfr1d,registerZ,registerI,WaistPos,InteWaist,registerARZI,register_fwhm_X,register_fwhm_Y,register_mesh]




def Source_slice(wfrG=[], zin = 0.0, zfin = -40., zstep=1.0):

    # ****************************** Finds the source position for slice slicen ************************************

    OUT    = findsource(wfrG, zin, zfin, zstep)
    registerZ       = OUT[1]
    registerI       = OUT[2]
    register_fwhm_X = OUT[6]
    register_fwhm_Y = OUT[7]
    writeres(res_dir+'IvsZ_sl'+str(slicen+1)+'.dat',registerZ,registerI)
    writeres(res_dir+'fwhmXvsZ_sl'+str(slicen+1)+'.dat',registerZ,register_fwhm_X)
    writeres(res_dir+'fwhmYvsZ_sl'+str(slicen+1)+'.dat',registerZ,register_fwhm_Y)

    WfrGTempl = OUT[0]

    arIP = array('f', [0]*WfrGTempl.mesh.nx*WfrGTempl.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIP, WfrGTempl, 6, 0, 3, WfrGTempl.mesh.eStart, 0, 0) #extracts intensity    
    arPP = array('d', [0]*WfrGTempl.mesh.nx*WfrGTempl.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
    srwl.CalcIntFromElecField(arPP, WfrGTempl, 0, 4, 3, WfrGTempl.mesh.eStart, 0, 0) 
    plotfield("Source_sl"+str(sample+1), WfrGTempl.mesh.xStart, WfrGTempl.mesh.xFin, WfrGTempl.mesh.nx, WfrGTempl.mesh.yStart, WfrGTempl.mesh.yFin, WfrGTempl.mesh.ny, arIP, arPP)

  


def xyprofave(pulse, tgrid, slicein, slicefin):
   
    print('Finding the average profile')
    for j in range(slicein,slicefin):
        slices = pulse[j]
        OUT = xyprofile(slices,tgrid)
        if j == slicein :
            Xcutave = OUT[0][1]
            Ycutave = OUT[1][1]
        else :
            Xcutave = Xcutave + OUT[0][1]
            Ycutave = Ycutave + OUT[1][1]
    Xcoord = OUT[0][0]
    Ycoord = OUT[1][0]
    Xcutave = Xcutave/(slicefin-slicein)
    Ycutave = Ycutave/(slicefin-slicein)

    return [[Xcoord, Xcutave],[Ycoord, Ycutave]]





def xydivave(pulse, tgrid, slicein, slicefin, lambd):
   
    print('Finding the average divergence')
    for j in range(slicein,slicefin):
        slices = pulse[j]
        OUT = xydiv(slices,tgrid,lambd)
        if j == slicein :
            Xcutave = OUT[0][1]
            Ycutave = OUT[1][1]
        else :
            Xcutave = Xcutave + OUT[0][1]
            Ycutave = Ycutave + OUT[1][1]
    Xcoord = OUT[0][0]
    Ycoord = OUT[1][0]
    Xcutave = Xcutave/(slicefin-slicein)
    Ycutave = Ycutave/(slicefin-slicein)

    return [[Xcoord, Xcutave],[Ycoord, Ycutave]]





def xyprofile(slices, tgrid):

    Xpos = int(np.argmax(np.abs(slices))/len(slices))
    Ypos = len(slices) -1 - int(np.argmax(np.abs(slices))/len(slices))
    ##print('max = ',np.max(np.abs(slices)),' maxpos = ',np.argmax(slices))
    ##print('XPOS = ',Xpos,' YPOS = ',Ypos,' MIDDLE = ',int(len(slices)/2))
    ##print('slmax = ',slices[Xpos,Ypos],' ',slices[Ypos,Xpos])
    #Xcut = np.abs(slices[int(len(slices)/2),:])**2
    #Ycut = np.abs(slices[:,int(len(slices)/2)])**2
    Xcut = np.abs(slices[Xpos,:])**2
    Ycut = np.abs(slices[:,Ypos])**2
    Xcoord = np.linspace(-tgrid,tgrid,len(Xcut))
    Ycoord = np.linspace(-tgrid,tgrid,len(Xcut))    
    #matplotlib.pyplot.figure()
    #matplotlib.pyplot.plot(Xcoord,Xcut)
    #matplotlib.pyplot.show()
    #matplotlib.pyplot.figure()
    #matplotlib.pyplot.plot(Ycoord,Ycut)
    #matplotlib.pyplot.show()

    return [[Xcoord, Xcut],[Ycoord, Ycut]]





def xydiv(slices, tgrid, lambd):    

    fslices = np.fft.fftshift(np.fft.fft2(slices))

    KXpos = int(np.argmax(np.abs(fslices))/len(fslices))
    KYpos = len(fslices) - 1 - int(np.argmax(np.abs(fslices))/len(fslices))
    ##print('KXPOS = ',KXpos,' KYPOS = ',KYpos,' MIDDLE = ',int(len(fslices)/2))
                  
    #KXcut = np.abs(fslices[int(len(slices)/2),:])**2
    #KYcut = np.abs(fslices[:,int(len(slices)/2)])**2
    KXcut = np.abs(fslices[KXpos,:])**2
    KYcut = np.abs(fslices[:,KYpos])**2  
    KXmax = 2*pi/(2*tgrid/len(KXcut))/2.0 #Total range is 2pi/(smallest size); smallest size is 2tgrid/len(KXcut); Max is total range/2
    KYmax = 2*pi/(2*tgrid/len(KYcut))/2.0
    THXmax = KXmax * lambd/(2*np.pi)
    THYmax = KYmax * lambd/(2*np.pi)        
    THXcoord = np.linspace(-THXmax,THXmax,len(KXcut))
    THYcoord = np.linspace(-THYmax,THYmax,len(KYcut))    
    #matplotlib.pyplot.figure()
    #matplotlib.pyplot.plot(THXcoord,KXcut)
    #matplotlib.pyplot.show()
    #matplotlib.pyplot.figure()
    #matplotlib.pyplot.plot(THYcoord,KYcut)
    #matplotlib.pyplot.show()
    
    return [[THXcoord, KXcut],[THYcoord, KYcut]]



def ProDiv_slice(sl = [], dgrid = 1.0, xlamds = 1.0):

    # ****************************** Calculates profile and divergence of slice slicen ************************************

    xy   = xyprofile(sl,dgrid)
    writeres(res_dir+'X_prof_sl'+str(slicen+1)+'.dat',xy[0][0],xy[0][1])
    writeres(res_dir+'Y_prof_sl'+str(slicen+1)+'.dat',xy[1][0],xy[1][1])
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(xy[0][0],xy[0][1])
    savefig(res_dir+'X_prof_sl'+str(slicen+1)+'.png')
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(xy[1][0],xy[1][1])
    savefig(res_dir+'Y_prof_sl'+str(slicen+1)+'.png')

    thxy = xydiv(sl,dgrid,xlamds)
    writeres(res_dir+'X_div_sl'+str(slicen+1)+'.dat',thxy[0][0],thxy[0][1])
    writeres(res_dir+'Y_div_sl'+str(slicen+1)+'.dat',thxy[1][0],thxy[1][1])
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(thxy[0][0],thxy[0][1])
    savefig(res_dir+'X_div_sl'+str(slicen+1)+'.png')
    matplotlib.pyplot.figure()
    matplotlib.pyplot.plot(thxy[1][0],thxy[1][1])
    savefig(res_dir+'Y_div_sl'+str(slicen+1)+'.png')

def GenWfr_slice(wfrG = []):

    # **********************Calculating Genesis Wavefront and extracting Intensity for slice slicen:

    XsG=wfrG.mesh.xStart
    XfG=wfrG.mesh.xFin
    nxG=wfrG.mesh.nx
    YsG=wfrG.mesh.yStart
    YfG=wfrG.mesh.yFin
    nyG=wfrG.mesh.ny

    arIG = array('f', [0]*nxG*nyG) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIG, wfrG, 6, 0, 3, wfrG.mesh.eStart, 0, 0) #extracts intensity
    arPG = array('d', [0]*nxG*nyG) #"flat" array to take 2D phase data (note it should be 'd')
    srwl.CalcIntFromElecField(arPG, wfrG, 0, 4, 3, wfrG.mesh.eStart, 0, 0) #extracts radiation phase
    plotfield("Gen_out_sl"+str(slicen+1), XsG, XfG, nxG, YsG, YfG, nyG, arIG, arPG)



def AveSourcePos(field = [], wfrG=[], zin = -20.0, zfin = -40.0, step = 1.0, slicein=-20.0, slicefin=-40.0):


    # ****************************** Finds the average source position ************************************

    OUTAVE = findsourceave(field, wfrG, -20.0, -40.0, 1.0, slicein, slicefin)
    registerZ_ave       = OUTAVE[0]
    registerI_ave       = OUTAVE[1]
    register_fwhm_X_ave = OUTAVE[4]
    register_fwhm_Y_ave = OUTAVE[5]
    reg_arzi_ave        = OUTAVE[6]
    writeres(res_dir+'IvsZ_ave.dat',registerZ_ave,registerI_ave)
    writeres(res_dir+'fwhmXvsZ_ave.dat',registerZ_ave,register_fwhm_X_ave)
    writeres(res_dir+'fwhmYvsZ_ave.dat',registerZ_ave,register_fwhm_Y_ave)
    maxind = np.argmax(registerI_ave)

    x = OUTAVE[7]
    y = OUTAVE[8]
    arIpy=np.sum(reg_arzi_ave[maxind],1)
    arIpx=np.sum(reg_arzi_ave[maxind],0)
    #These are summed! The previous ones are cut 
    writeres(res_dir+'IvsX_ave.dat',x,arIpx)
    writeres(res_dir+'IvsY_ave.dat',y,arIpy)

    return OUTAVE[2]   


      
def prtfluence(wfr = [], arI=[], MESH = [] ):

    Xs = MESH[0]
    Xf = MESH[1]
    nx = MESH[2]
    Ys = MESH[3]
    Yf = MESH[4]        
    ny = MESH[5]
    pixel_area=(Xf-Xs)*(Yf-Ys)/nx/ny
    integral=sum([x*pixel_area for x in arI])
    np.disp('integral_area(a.u.)= '+str(integral/1000000))

    np.disp('max_intensity= '+str(round_sig(max(arI),3)))
    np.disp('pixel_area= '+str(pixel_area))
    Ntotph = 1
    flux = fluence(wfr.arEx, Xf-Xs, Yf-Ys, nx, ny, 1)
    print('X size = ',Xf-Xs,' m')
    print('Y size = ',Yf-Ys,' m')
    print('X mesh = ',nx,' points')
    print('Y mesh = ',ny,' points')
    print('Total fluence is ',flux,'*Nph in ph/m2')
    print('also corresponding to ',flux/1.0e6,'*Nph in ph/mm2')


def AVE_prtfluence(INT_arI_ave_P = [], MESH = []):

    XsP0 = MESH[0]
    XfP0 = MESH[1]
    nxP0 = MESH[2]
    YsP0 = MESH[3]
    YfP0 = MESH[4]        
    nyP0 = MESH[5]
    pixel_area=(XfP0-XsP0)*(YfP0-YsP0)/nxP0/nyP0
    integral=sum([x*pixel_area for x in INT_arI_ave_P])    
    np.disp('integral_area(a.u.)= '+str(integral/1000000))

    np.disp('max_intensity= '+str(round_sig(max(INT_arI_ave_P),3)))
    np.disp('pixel_area= '+str(pixel_area))
    Ntotph = 1
    flux = fluenceI(INT_arI_ave_P, XfP0-XsP0, YfP0-YsP0, nxP0, nyP0, 1)
    print('X size = ',XfP0-XsP0,' m')
    print('Y size = ',YfP0-YsP0,' m')
    print('X mesh = ',nxP0,' points')
    print('Y mesh = ',nyP0,' points')
    print('Total fluence is ',flux,' ph/m2')
    print('also corresponding to ',flux/1.0e6,' ph/mm2')


'''

'''
#################################################################################################
####################           Routines strictly used in prop_par.py        #####################
#################################################################################################

def GeneRead(comm=[], filename='D:\prog\data\PROP\sase2.4p1keV.tap.1.out.dfl', Tmesh=301, PhEne = 4100.0, dgrid = 1.0):
    """
    Reads the Genesis file and gives back wavefront backbone and field as a np.array
    """
    
    
    wfrG = SRWLWfr()               #Initial Electric Field Wavefront
    wfrG.allocate(1, Tmesh, Tmesh) #Numbers of points vs Photon Energy (1), Horizontal and Vertical Positions (dummy)
    wfrG.mesh.zStart =  0.0        #Longitudinal Position [m] at which Electric Field has to be calculated, i.e. the position of the first optical element
    wfrG.mesh.eStart =  PhEne      #Initial Photon Energy [eV]
    wfrG.mesh.eFin   =  PhEne      #Final Photon Energy [eV]
    wfrG.mesh.xStart = -dgrid      #Initial Horizontal Position [m]
    wfrG.mesh.xFin   =  dgrid      #Final Horizontal Position [m]
    wfrG.mesh.yStart = -dgrid      #Initial Vertical Position [m]
    wfrG.mesh.yFin   =  dgrid      #Final Vertical Position [m]

    wfrG.partBeam.partStatMom1.x  = 0.0     #Some information about the source in the Wavefront structure
    wfrG.partBeam.partStatMom1.y  = 0.0
    wfrG.partBeam.partStatMom1.z  = 0.0
    wfrG.partBeam.partStatMom1.xp = 0.0
    wfrG.partBeam.partStatMom1.yp = 0.0

    #Reading Genesis field file
    print('Reading Genesis field file...')
    field=readRadiationFile_mpi(comm=comm, fileName=filename, npoints=Tmesh)
    print('...Done.',rank,'. ',len(field))
    
    return [wfrG, field]


def Source_Prop(Len_drft = 0.1, Wfrin = []):
    
    #**********************Calculating virtual Source Wavefront and extracting Intensity:

    Drift_1 = SRWLOptD(Len_drft) #Drift    
    # 2colors U2 sase3 :
    #                    [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11] 
    ppDrift_1 =          [ 0,  0, 1.0,  0,  0, 1.0, 2.0, 1.0, 2.0,  0,  0,   0]   
    ppFinal =            [ 0,  0, 1.0,  0,  0, 0.5, 1.0, 0.5, 1.0,  0,  0,   0]

    # Usual source sase1
    #ppDrift_1 =          [ 0,  0, 1.0,  0,  0, 1.0, 2.0, 1.0, 2.0,  0,  0,   0]   
    #ppFinal =            [ 0,  0, 1.0,  0,  0, 0.25, 2.0, 0.25, 2.0,  0,  0,   0]
   
    
    #2colors_backUndu sase3
    #ppDrift_1 =          [ 0,  0, 1.0,  0,  0, 16.0, 1.0, 16.0, 1.0,  0,  0,   0]   
    #ppFinal =            [ 0,  0, 1.0,  0,  0, 1.0, 0.1, 1.0, 0.1,  0,  0,   0]
    
    
    optBL = SRWLOptC([Drift_1],  [ppDrift_1, ppFinal])
    wfr1d = copy.deepcopy(Wfrin)
    srwl.PropagElecField(wfr1d, optBL)
    arI = array('f', [0]*wfr1d.mesh.nx*wfr1d.mesh.ny) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arI, wfr1d, 6, 0, 3, wfr1d.mesh.eStart, 0, 0) #extracts intensity    
    arP = array('d', [0]*wfr1d.mesh.nx*wfr1d.mesh.ny) #"flat" array to take 2D phase data (note it should be 'd')
    srwl.CalcIntFromElecField(arP, wfr1d, 0, 4, 3, wfr1d.mesh.eStart, 0, 0) #extracts radiation phase

    return [wfr1d, arI, arP]
    
    

def FullProp(zprop = 933.0, wfrS = [], line = SASE1_SPB_line, ppFinal = [ 0,  0, 1.0,  1,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]):
    
    print line
    [ Act_Elem_box, Act_Param_box ] = line(zprop)
    Act_Param_box.append( ppFinal )
    OptBL = SRWLOptC(Act_Elem_box, Act_Param_box)
    wfrP = copy.deepcopy(wfrS) #Genesis backpropagated    
    srwl.PropagElecField(wfrP, OptBL)
    
    return wfrP
    
    
def AVE_Save_Inte(wfrP=[], MESH = [], INT_arI_ave_P = [],  firstslice = 1):

    XsP=wfrP.mesh.xStart
    XfP=wfrP.mesh.xFin
    nxP=wfrP.mesh.nx
    YsP=wfrP.mesh.yStart
    YfP=wfrP.mesh.yFin
    nyP=wfrP.mesh.ny
    
    arIP = array('f', [0]*nxP*nyP) #"flat" array to take 2D intensity data
    srwl.CalcIntFromElecField(arIP, wfrP, 6, 0, 3, wfrP.mesh.eStart, 0, 0) #extracts intensity
    if firstslice == 1:        
        INT_arI_ave_P = np.array(arIP)        
        XsP0 = XsP
        XfP0 = XfP
        YsP0 = YsP
        YfP0 = YfP
        nxP0 = nxP
        nyP0 = nyP
    else:
        XsP0 = MESH[0]
        XfP0 = MESH[1]
        nxP0 = MESH[2]
        YsP0 = MESH[3]
        YfP0 = MESH[4]        
        nyP0 = MESH[5]
        XMESH = np.linspace(XsP0, XfP0, nxP0)
        YMESH = np.linspace(YsP0, YfP0, nyP0)
        dataI = np.array( arIP )
        zI    = dataI.reshape( nyP, nxP )
        Xact = np.linspace(XsP, XfP, nxP)
        Yact = np.linspace(YsP, YfP, nyP)
        Interdata = interpolate.interp2d(Xact, Yact, zI)
        zI_inte    = Interdata(XMESH, YMESH)
        INT_arIP   = zI_inte.reshape(nyP*nxP)
        INT_arI_ave_P = INT_arI_ave_P + np.array(INT_arIP)

    wfrP.mesh.xStart = XsP0
    wfrP.mesh.xFin   = XfP0
    wfrP.mesh.nx     = nxP0
    wfrP.mesh.yStart = YsP0
    wfrP.mesh.yFin   = YfP0
    wfrP.mesh.ny     = nyP0
    wfrP.mesh.ne     = 1
    

    return [INT_arI_ave_P, [XsP0,XfP0,nxP0,YsP0,YfP0,nyP0]]


######################################
####GRATINGS ROUTINES TO BE TESTED####
######################################


def fftpulse(field_in = []):

    fft_field = field_in
    nslice = len(field_in)
    xmesh  = len(field_in[0])
    ymesh  = len(field_in[0,0])

    for j in range(xmesh):
        for k in range(ymesh):
            fft_field[:,j,k] = np.fft.fft(field_in[:,j,k])    

    return fft_field

def ifftpulse(field_in = []):

    ifft_field = field_in
    nslice = len(field_in)
    xmesh  = len(field_in[0])
    ymesh  = len(field_in[0,0])

    for j in range(xmesh):
        for k in range(ymesh):
            ifft_field[:,j,k] = np.fft.ifft(field_in[:,j,k])    


    return ifft_field

    

def grating(field_in = [], leng = 1, xlamds = 0.0, z1 = 0.0, z2 = 0.0, r_t = 0.0, r_s = 0.0, theta_i = 0.0, theta_d =0.0, dtheta_d, d0 = 0.0, d1 = 0.0, d2 = 0.0, aberration = False):
    """
    % field_in is the input complex field
    % field_out is the output complex field
    % leng is transverse size of the mesh (dgrid*2)
    % xlamds is wavelength
    % z1 and z2 are distances to source and image in order to calculate 
    %   aberration effect (we use semi-analytical approach to propagate, 
    %   so we need to provide this information)
    % R_t is tangential radius of curvature
    % R_s is sagittal radius of curvature
    % Teta_i and Teta_d are incidence and reflection angles (obviously one can 
    %   calculate it now internaly within function. I keep dummy input variable
    %   to save time on revriting all versiong of script code.
    % dTheta_d is angle deviation of the field from principle ray direction
    %   after grating. It is used to introduce different angles of reflected
    %   field from grating at different  wavelengths if time-dependent (3D) field
    %   is propagated, resulting in pusle-front tilt.
    % D0, D1, D2 are grating coefficients (in out case 1123[l/mm] 1.6[l/mm2] 0.002[l/mm3])
    % aberration is logical value, that defines whether include aberration effects, or not
    """

    
    theta_d  = np.arccos(np.cos(theta_i)-xlamds*d0*1.0e3)
    k        = 2.0*np.pi/xlamds
    gr_asym  = theta_i/theta_d         # grating asymmetry factor GG: assumes grazing angles sin angle ~ angle (?)
    gr_sigma = 1.0/(d0*1.0e3)          # sigma parameter of grating (1/k)
    gr_alpha = d1*1.0e6/(d0*1.0e3)**2  # alpha parameter of grating (old=1,8) (n1/k/k) !!
    d2       = d2*1.0e9                # in 1/m**3
    if d2 == 0.0:
        n=0
        d2=1.0
    else
        n=1

    # calculation of the sagittal focal distance 
    f_s = r_s/( np.sin(theta_i) + np.sin(theta_d) )
    #print 'f_g_s = ', f_s

    # calculation of the tangential focal distance 
    #1. vls contribution
    f_t_vls   = gr_sigma**2 * np.sin(theta_d)**2/xlamds/gr_alpha  #!!! original!!!!
    f_g_t_vls = gr_sigma**2/xlamds/gr_alpha
    #2. curvature contribution
    f_t_c     = r_t/(1.0/theta_d+theta_i/theta_d**2)              #tangential focal distance due to curvature
    #3. final tangential focal distance
    f_t       = 1.0 / ( 1.0/f_t_vls + 1.0/f_t_c )    
    #print 'f_g_t = ',f_t,' m,  f_g_t_c = ',f_t_c,' m,  f_g_t_vls = ',f_t_vls,' m'

    m      = len(fiel_in)
    dx     = leng/m
    dy     = dx
    x1     = np.linspace( (m-1.0)/2.0*dx, ((m-1.0)/2.0 + 1-m))*dx, m )
    y1     = np.linspace( (m-1.0)/2.0*dy, ((m-1.0)/2.0 + 1-m))*dy, m )
    xx, yy = np.meshgrid(x1, y1) 

    #c_20=(theta_i**2/z1 + theta_i/r_t + theta_d**2/z2 + theta_d/r_t)/2
    if aberration == True:
        c_30 = ( n*np.pi*2.0/K/3.0*d2 ) + ( np.sin(theta_i)**2/z1 - np.sin(theta_i)/r_t ) * np.cos(theta_i)/2.0/z1 - ( np.sin(theta_d)**2/z2 - sin(theta_d)/r_t )*cos(theta_d)/2.0/z2
        c_12 = ( -theta_i/r_s/z1 + 1.0/z1**2 + theta_d/r_s/z2 - 1/z2**2)/2.0
    else
        c_30 = 0
        c_12 = 0
     

    # enegry conservation calculations
    dm   = np.round( m*(1.0-gr_asym)/2.0 )
    p_in = np.sum(np.abs(field_in(dm:m-dm))**2)
    
    fx    = interpolate.interp2d(xx, yy, field_in)
    x     = fx(xx, yy*gr_asym)
    p_out = np.sum(np.abs(x)**2)
    
    if pout ~= 0.0:
        x = x*np.sqrt(p_in/p_out)
    
    x = x * np.exp(-1j*k/f_s/2.0*(xx**2))                #focusing sag
    x = x * np.exp(-1j*k/f_t/2.0*(yy**2))                #focusing tang
    x = x * np.exp(1j*k*c_30*(yy**33)/theta_d**3)        #aberration c_30 
    x = x.* np.exp(1j*k*c_12*xx**2 * yy/theta_d)         #aberration c_12

    x = np.abs(x)*np.exp(1j*(np.angle(x)+dtheta_d*k*yy)) #phase tilt
        
    field_out=x;
        
    return field_out
'''
