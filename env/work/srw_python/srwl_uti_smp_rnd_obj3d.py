#############################################################################
# uti_sim module: misc. sample simulation utilities / functions
# v 0.01
# Authors: Himanshu Goel (SBU/ECE), O.C. (BNL/NSLS-II)
#############################################################################

from __future__ import print_function #Python 2.7 compatibility
from array import *
from math import *
from copy import *
import random

# ********************** Simulate Brownian motion of 3D (Nano-) Objects
def setup_list_obj3d(_n, _ranges, _cen, _dist, _obj_shape, _allow_overlap=True, _seed=None, _fp=None, _sep='\t', _nl='\n'): #OC26082021
#def setup_list_obj3d(_n, _ranges, _cen, _dist, _obj_shape, _allow_overlap=True, _seed=0, _fp=None, _sep='\t', _nl='\n'): #OC16082021
#def setup_sph3d_from_distribution(nElements, szDist, posDist, volume, allowOverlap=False, seed=0, out_file=None, delim='\t', newline='\n'):
    """
    Creates a list of 3D objects (randomly) distributed over a volume.
    :param _n: number of 3D objects to generate;
    :param _ranges: list of ranges of positions x, y, z (in [m]), defining the volume within which the 3D objects should be generated;
    :param _cen: list of Cartesian coordinates x, y, z (in [m]) of the center point of the volume within which the 3D objects objects should be generated;
    :param _dist: (list of) parameters defining object distribution in 3D space; the following are the supported distributions and parameters:
            'uniform' or ['uniform'] for the uniform distribution of objects within the the volume defined by _dims and _cen;
            ['gauss', sigx, sigy, sigz, xc, yc, zc] for the Gaussian distribution with the standard deviations vs x, y, z defined by sigx, sigy, sigz and center position defined by xc, yc, zc (if the latter values are not supplied the center of the distribution will be defined by _cen);
    :param _obj_shape: list of 3D object shape parameters, with its elements defining:
            [0]: object type (currently only 'S'- sphere is supported);
            [1]: type of object size distribution (currently 'uniform', 'gauss' and 'exp' are supported); the following elements contain parameters of the object size distribution;
            - in the case of 'uniform' distribution:
            [2]: is the minimum size of 3D object (diameter for sphere);
            [3]: is the maximum size of 3D object (diameter for sphere);
            - in the case of Gaussian distribution ('gauss'):
            [2]: is the standard deviation of the 3D object size (diameter for sphere);
            - in the case of exponential distribution ('exp'):
            [2]: is the mean value of the 3D object size (diameter for sphere);
            subsequent parameters may be used to define orientation of non-uniform 3D objects in space;
            when object shapes other than spheres will be supported, _obj_shape may be defined as 2D list, with each 1D list defining parameters of one type of 3D objects;
    :param _allow_overlap: whether or not to allow the 3D object to overlap;
    :param _seed: random number generator seed value;
    :param _out_file: optional name of the file to save the generated list to (if not provided, the generated list won't be saved to a file);
    :param _delim: delimiter between columns of the file;
    :param _newline: delimiter between rows of the file;
    :return: 2D list containing 3D object definitions.
    """
    
    def sample_dist_from_name(_par): #OC16082021 (modified, the original version was by #HG)
        """
        Samples a distribution given a random.Random instance and distribution name + parameters.
        :param rng: random.Random instance.
        :param params: list of distribution name + parameters.
        :return: Sampled value.
        """
        if type(_par[0]) == float: #Assume a constant value if no distribution is specified
            return _par[0]
        if((_par[0] == 'gauss') or (_par[0] == 'normal')):
            return random.gauss(_par[1], _par[2])
        #elif((_par[0] == 'exponential') or (_par[0] == 'exp')):
        #    return random.expovariate(_par[1])
        elif _par[0] == 'uniform':
            return random.uniform(_par[1], _par[2])
        return None

    if(_seed is not None): random.seed(_seed) #OC27082021
    #random.seed(_seed)

    min_pos = [_cen[j] - 0.5*_ranges[j] for j in range(3)]
    max_pos = [_cen[j] + 0.5*_ranges[j] for j in range(3)]

    posDistPar = None
    if(isinstance(_dist, list)): 
        if(_dist[0] == 'uniform'): posDistPar = [[_dist[0], min_pos[j], max_pos[j]] for j in range(3)]
        elif((_dist[0] == 'gauss') or (_dist[0] == 'normal')): 
            if(len(_dist) == 4): posDistPar = [[_dist[0], _cen[j], _dist[j + 1]] for j in range(3)]
            elif(len(_dist) == 7): posDistPar = [[_dist[0], _dist[j + 4], _dist[j + 1]] for j in range(3)]
    else: posDistPar = [[_dist, min_pos[j], max_pos[j]] for j in range(3)]

    if(posDistPar is None): raise Exception("Incorrect definition of 3D object positions")

    nObjShapes = 1
    if(isinstance(_obj_shape[0], list)): 
        nObjShapes = len(_obj_shape)

    curObjType = None
    curObjSizeDistPar = None

    if(nObjShapes <= 1):
        curObjType = _obj_shape[0]
    
    if(curObjType is not None):
        if(curObjType == 'S'):
            if((_obj_shape[1] == 'uniform') or (_obj_shape[1] == 'gauss') or (_obj_shape[1] == 'normal')): #OC26082021
                curObjSizeDistPar = [_obj_shape[1], _obj_shape[2], _obj_shape[3]]
            #if(_obj_shape[1] == 'uniform'): 
            #    curObjSizeDistPar = [_obj_shape[1], _obj_shape[2], _obj_shape[3]]
            #elif((_obj_shape[1] == 'gauss') or (_obj_shape[1] == 'normal')):
            #    curObjSizeDistPar = [_obj_shape[1], 0., _obj_shape[2]]
        #elif(curObjType == ...):

    of = None
    if _fp is not None: of = open(_fp, 'w')

    resList = []
    for i in range(_n):

        if(nObjShapes > 1):
            curObjShape = _obj_shape[random.randint(0, nObjShapes - 1)]
            curObjType = curObjShape[0]
            if(curObjType == 'S'):
                if((curObjShape[1] == 'uniform') or (curObjShape[1] == 'gauss') or (curObjShape[1] == 'normal')): 
                    curObjSizeDistPar = [curObjShape[1], curObjShape[2], curObjShape[3]]
                #if(curObjShape[1] == 'uniform'): 
                #    curObjSizeDistPar = [curObjShape[1], curObjShape[2], curObjShape[3]]
                #elif((curObjShape[1] == 'gauss') or (curObjShape[1] == 'normal')):
                #    curObjSizeDistPar = [curObjShape[1], 0., curObjShape[2]]
            #elif(curObjType == ...):

        curObjDescr = None
        while True:
            curObjDescr = [sample_dist_from_name(posDistPar[j]) for j in range(3)] #Coordinates of the current object (candidate)

            if(curObjType == 'S'):
                while True: #OC26082021 (?)
                    curObjD = sample_dist_from_name(curObjSizeDistPar)
                    if(curObjD > 0): break
                curObjR = 0.5*curObjD

                #Checking if the object crosses boundary of the volume:
                if all([(curObjDescr[j] >= min_pos[j]) and (curObjDescr[j] < max_pos[j]) for j in range(3)]): #OC26082021 (to be able to create consistent "slices" for "thick" samples)
                #if all([(curObjDescr[j] - curObjR >= min_pos[j]) and (curObjDescr[j] + curObjR <= max_pos[j]) for j in range(3)]): 

                    overlapFound = False
                    if(not _allow_overlap): #Checking if the object overlaps with any other object in the volume:
                        if(curObjType == 'S'):
                            for ii in range(len(resList)):
                                exObjDescr = resList[ii]
                                if(exObjDescr[3] == 'S'):
                                    dx = curObjDescr[0] - exObjDescr[0]
                                    dy = curObjDescr[1] - exObjDescr[1]
                                    dz = curObjDescr[2] - exObjDescr[2]
                                    twoR = curObjR + exObjDescr[4]
                                    if(dx*dx + dy*dy + dz*dz < twoR*twoR):
                                        overlapFound = True
                                        break

                    if(not overlapFound):
                        curObjDescr.append(curObjType)

                        if(curObjType == 'S'):
                            curObjDescr.append(curObjR) #To check if this should be radius or diameter
                            resList.append(curObjDescr)
                            if of is not None:
                                of.write('{:15.11g}{:}{:15.11g}{:}{:15.11g}{:}{:}{:}{:15.11g}{:}'.format(curObjDescr[0], _sep, curObjDescr[1], _sep, curObjDescr[2], _sep, curObjDescr[3], _sep, curObjDescr[4], _nl))
                                #of.write('{:15.12e}{:}{:15.12e}{:}{:15.12e}{:}S{:}{:15.12e}{:}'.format(curObjDescr[0], _sep, curObjDescr[1], _sep, curObjDescr[2], _sep, curObjR, _nl))
                        break

    if of is not None: of.close()

    #rng = random.Random(seed)

    # # assume volume is centered at 0,0,0
    # min_pos = [volume[x] * -0.5 for x in range(3)]
    # max_pos = [volume[x] * 0.5 for x in range(3)]

    # # if position distribution info is single length, the same distribution parameters are used for all three positions
    # if len(posDist) == 1:
    #     posDist = [posDist[0] for _ in range(3)]

    # # don't write to file if a name hasn't been specified.
    # if out_file != None:
    #     of = open(out_file, 'w')

    # defns = []
    # for _ in range(nElements):
    #     while True:
    #         r = abs(sample_dist_from_name(rng, szDist))
    #         xyz = [sample_dist_from_name(rng, posDist[x]) for x in range(3)]
    #         if all([xyz[x] - r > min_pos[x] and xyz[x] + r < max_pos[x] for x in range(3)]):
    #             break  # Draw samples until the object does not intersect with the volume

    #         if not allowOverlap:  # check for overlap (inefficiently), to be more efficient and accurate the C++ implementation can use an octree, or just a regular grid with bins
    #             overlapDetected = False
    #             for d_xyz0, d_xyz1, d_xyz2, _, d_r in defns:
    #                 dist = (d_xyz0 - xyz[0]) ** 2 + (d_xyz1 - xyz[1]) ** 2 + (d_xyz2 - xyz[2]) ** 2 # compute distance squared to the current sphere
    #                 if dist < (d_r + r) ** 2:
    #                     overlapDetected = True
    #                     break
    #             if not overlapDetected: # no overlap found, this sample can be accepted
    #                 break

    #     defns.append([xyz[0], xyz[1], xyz[2], 'S', r])

    #     if out_file != None:
    #         of.write('{:15.12f}{:}{:15.12f}{:}{:15.12f}{:}S{:}{:15.12f}{:}'.format(xyz[0], delim, xyz[1], delim, xyz[2], delim, delim, r, newline))

    # if out_file != None:
    #     of.close()

    return resList

# ********************** Simulate Brownian motion of 3D (Nano-) Objects
def brownian_motion3d(_obj_crd, _viscosity=1.e-3, _temperature=293, _timestep=0.1, _duration=10, _seed=None, _fp=None, _sep='\t', _nl='\n'): #OC08102021
#def brownian_sim3d(shapeData, _viscosity = 1e-3, _temperature = 293, _timestep = 0.1, _duration = 10, seed = None, out_dir=None, delim='\t', newline='\n'): #HG06302021
    """
    Apply 3D Brownian motion simulation to a set of input objects.
    :param _obj_crd: 2D list of initial coordinates of objects for which the motion is to be simulated; _obj_crd[i] should contain a list of Cartesian coordinates (x, y, z) of objects, possibly followed by additional information on the objects (that won't be modified)
    :param _viscosity: viscosity of the simulated medium [Pa*s]
    :param _temperature: temperature of the medium [K]
    :param _timestep: simulation timestep [s]
    :param _duration: length of the simulation [s]
    :return: list of shapes and their positions throughout the simulation
    """

    try:
        import numpy as np
    except:
        raise Exception('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following command to install it: \npip install numpy')

    k_B = 1.380649e-23 #J/K
    N = len(_obj_crd) #OC10082021
    #N = len(shapeData)
    _timestep_count = int(_duration / _timestep)

    D = [k_B * _temperature / (3 * np.pi * _viscosity * 2 * x[4]) for x in _obj_crd] #OC10082021
    #D = [k_B * _temperature / (3 * np.pi * _viscosity * 2 * x[4]) for x in shapeData]
    k = [sqrt(2 * D_v * _timestep) for D_v in D]

    rng = np.random.seed(seed = _seed) #OC31082021 (changing because of a problem running on Linux)
    dx = np.random.normal(scale = k, size = (_timestep_count, N))
    dy = np.random.normal(scale = k, size = (_timestep_count, N))
    dz = np.random.normal(scale = k, size = (_timestep_count, N))

    #rng = np.random.default_rng(seed = _seed)
    #dx = rng.normal(scale = k, size = (_timestep_count, N))
    #dy = rng.normal(scale = k, size = (_timestep_count, N))
    #dz = rng.normal(scale = k, size = (_timestep_count, N))

    poses = [_obj_crd] #OC10082021
    #poses = [shapeData]

    for step in range(1, _timestep_count + 1):
        cur_poses = []
        
        of = None if _fp is None else open(_fp%(step), 'w') #OC08102021
        #if out_dir != None:
        #    of = open(out_dir + '/%d.def'%(step), 'w')

        for i in range(N):

            #OC10082021: To keep other information about objects in that list without parsing it explicitly
            newLineObjCrd = copy(_obj_crd[i]) 
            newLineObjCrd[0] = poses[-1][i][0] + dx[step - 1][i]
            newLineObjCrd[1] = poses[-1][i][1] + dy[step - 1][i]
            newLineObjCrd[2] = poses[-1][i][2] + dz[step - 1][i]
            cur_poses.append(newLineObjCrd)

            #cur_poses.append([poses[-1][i][0] + dx[step - 1][i], 
            #                  poses[-1][i][1] + dy[step - 1][i], 
            #                  poses[-1][i][2] + dz[step - 1][i], 
            #                  shapeData[i][3], 
            #                  shapeData[i][4]])

            if of is not None: #OC17082021
                of.write('{:15.11g}{:}{:15.11g}{:}{:15.11g}{:}{:}{:}{:15.11g}{:}'.format(cur_poses[-1][0], _sep, cur_poses[-1][1], _sep, cur_poses[-1][2], _sep, cur_poses[-1][3], _sep, cur_poses[-1][4], _nl))
            #if out_dir != None:
            #    of.write('{:15.12f}{:}{:15.12f}{:}{:15.12f}{:}{:}{:}{:15.12f}{:}'.format(cur_poses[-1][0], delim, cur_poses[-1][1], delim, cur_poses[-1][2], delim, cur_poses[-1][3], delim, cur_poses[-1][4], newline))

        if of is not None: of.close() #OC17082021
        #if out_dir != None:
        #    of.close()

        poses.append(cur_poses)
    return poses
