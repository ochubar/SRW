# -*- coding: utf-8 -*-
#############################################################################
# SRW Samples library
# Authors/Contributors: Maksim Rakitin, Rebecca Ann Coles
# October 26, 2016
# April 11, 2020
#############################################################################

import math
import os
import sys
import time
from array import array

import srwlib
import uti_io


# ********************** The class for Samples:
class SRWLUtiSmp:
    """The class for Samples from image file (e.g., .tif), NumPy array (.npy), etc.

    :param file_path: full path to the image or the saved NumPy array.
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_show_images: a flag to show the initial and processed images.
    :param is_save_images: a flag to save the initial and processed images.
    :param raw_image_name: the name of the raw file in case if it's saved.
    :param processed_image_name: the name of the processed file in case if it's saved.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    """

    def __init__(self, file_path,
                 area=None, rotate_angle=0, rotate_reshape=True, cutoff_background_noise=0.25, background_color=0,
                 tile=None, shift_x=None, shift_y=None, invert=None,
                 is_show_images=False, is_save_images=False,
                 raw_image_name='raw', processed_image_name='processed', prefix='', output_image_format=None, #):
                 cutoff_max_fact=None, max_color=255): #OC29052020
        # Input parameters:
        self.file_path = file_path
        self.area = area
        self.rotate_angle = rotate_angle if rotate_angle is not None else 0
        self.rotate_reshape = True if rotate_reshape else False
        self.cutoff_background_noise = cutoff_background_noise if cutoff_background_noise is not None else 0
        self.background_color = background_color if background_color is not None else 0
        self.invert = invert
        self.tile = tile
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.is_show_images = is_show_images
        self.is_save_images = is_save_images
        output_image_format = os.path.splitext(file_path)[1].replace('.', '') if not output_image_format else \
            output_image_format
        self.raw_image_name = self._add_prefix(prefix, raw_image_name, output_image_format)
        self.processed_image_name = self._add_prefix(prefix, processed_image_name, output_image_format)

        self.cutoff_max_fact = cutoff_max_fact #OC29052020
        self.max_color = max_color if max_color is not None else 255 #OC29052020

        # Output parameters:
        self.data = None
        self.raw_image = None
        self.processed_image = None
        self.nx = None
        self.ny = None
        self.limit_value = None

        # Check input type automatically:
        self.input_type = self._check_input_type(file_path)

        # Set the dir where to save the files:
        self.save_dir = os.path.abspath(os.path.dirname(file_path))

        # Check the input file(s):
        self._check_files()

        # Process input:
        self.read_sample()

        # Show the resulted images:
        if self.is_show_images:
            self.show_images()

        # Save the resulted images:
        if self.is_save_images:
            self.save_images()

    def get_data_from_image(self):
        import_err_msg = '"{}" library cannot be imported. Please install it first with "pip install {}".'
        try:
            import numpy as np
        except ImportError:
            raise ImportError(import_err_msg.format('NumPy', 'numpy'))
        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format('PIL', 'pillow'))
        try:
            from scipy.ndimage.interpolation import rotate
        except ImportError:
            raise ImportError(import_err_msg.format('SciPy', 'scipy'))

        d = uti_io.read_image(image_path=self.file_path)
        self.data = d['data']
        self.raw_image = d['raw_image']
        self.limit_value = d['limit_value']

        # Remove background noise:
        #OC29052020 (moved down)
        #assert 0 <= self.background_color <= 255, 'Background color ({}) should be between 0 and 255'.format(
        #    self.background_color)
        #assert 0 <= self.cutoff_background_noise <= 1, 'Cutoff background noise ({}) should be between 0 and 1'.format(
        #    self.cutoff_background_noise)
        #self.data[np.where(self.data < self.limit_value * self.cutoff_background_noise)] = np.uint16(
        #    self.background_color)

        if self.area:
            assert type(self.area) in [tuple, list] and len(self.area) == 4, \
                'The area should be a list/tuple and contain 4 elements (x_start, x_end, y_start, y_end)'
            assert self.area[0] <= self.data.shape[1] and self.area[1] <= self.data.shape[1] and \
                   self.area[2] <= self.data.shape[0] and self.area[3] <= self.data.shape[0], \
                '\n  x_start ({}) and x_end ({}) should be less than {}\n  y_start ({}) and y_end ({}) should be less than {}'.format(
                    self.area[0], self.area[1], self.data.shape[1],
                    self.area[2], self.area[3], self.data.shape[0],
                )
            assert self.area[0] < self.area[1] and \
                   self.area[2] < self.area[3], \
                '\n  x_start ({}) should be less than x_end ({})\n  y_start ({}) should be less than y_end ({})'.format(
                    self.area[0], self.area[1],
                    self.area[2], self.area[3],
                )
            self.data = self.data[self.area[2]:self.area[3], self.area[0]:self.area[1]]

        # Remove background noise:
        #OC29052020 (moved here to make operations on reduced amount of data)
        assert 0 <= self.background_color <= 255, 'Background color ({}) should be between 0 and 255'.format(
            self.background_color)
        assert 0 <= self.cutoff_background_noise <= 1, 'Cutoff background noise ({}) should be between 0 and 1'.format(
            self.cutoff_background_noise)
        self.data[np.where(self.data < self.limit_value * self.cutoff_background_noise)] = np.uint16(
            self.background_color)

        if(hasattr(self, 'cutoff_max_fact')): #OC29052020
           if(self.cutoff_max_fact is not None):
                assert 0 <= self.cutoff_max_fact <= 1, 'Maximal cutoff signal value ({}) should be between 0 and 1'.format(self.cutoff_max_fact)
                assert 0 <= self.max_color <= 255, 'Maximum color value ({}) should be between 0 and 255'.format(self.max_color)
                self.data[np.where(self.data >= self.limit_value * self.cutoff_max_fact)] = np.uint16(self.max_color)

        if self.tile:
            assert type(self.tile) in [list, tuple], 'The type of tile ({}) should be list/tuple'.format(
                type(self.tile).__name__)
            assert len(self.tile) == 2, 'The size ({}) of the list/tuple should be 2'.format(len(self.tile))
            self.data = np.tile(self.data, self.tile)

        if self.rotate_angle:
            assert -360 < self.rotate_angle < 360, 'The angle should be from -360 to 360 degrees.'
            self.data = rotate(self.data, self.rotate_angle, axes=(1, 0), reshape=self.rotate_reshape, output=None,
                               order=0, mode='constant', cval=self.background_color, prefilter=False)

        if self.shift_x:
            assert type(self.shift_x) is int, 'Type of shift_x ({}) should be int'.format(type(self.shift_x).__name__)
            if self.shift_x > 0:
                self.data = np.pad(self.data, ((0, 0), (self.shift_x, 0)), mode='constant',
                                   constant_values=(self.background_color))[:, :-self.shift_x]
            else:
                self.data = np.pad(self.data, ((0, 0), (0, -self.shift_x)), mode='constant',
                                   constant_values=(self.background_color))[:, -self.shift_x:]

        if self.shift_y:
            assert type(self.shift_y) is int, 'Type of shift_y ({}) should be int'.format(type(self.shift_y).__name__)
            if self.shift_y < 0:
                self.data = np.pad(self.data, ((-self.shift_y, 0), (0, 0)), mode='constant',
                                   constant_values=(self.background_color))[:self.shift_y, :]
            else:
                self.data = np.pad(self.data, ((0, self.shift_y), (0, 0)), mode='constant',
                                   constant_values=(self.background_color))[self.shift_y:, :]

        if self.invert:
            self.data = np.invert(self.data)

        self.processed_image = Image.fromarray(self.data)

    def get_image_from_data(self):
        data = np.load(self.file_path)
        self.limit_value = 255
        #Suggested by RAC(?):
        self.data = (np.subtract(data, data.min(), dtype=np.float32)) / (np.subtract(data.max(), data.min(), dtype=np.float32)) * self.limit_value
        #self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        self.raw_image = Image.fromarray(self.data.astype(np.uint8))
        self.processed_image = self.raw_image

    def read_sample(self):
        if self.input_type == 'image':
            self.get_data_from_image()
        elif self.input_type == 'npy':
            self.get_image_from_data()
        else:
            raise NotImplementedError(
                'Processing of the "{}" input type is not implemented yet.'.format(self.input_type))
        self.nx = self.data.shape[1]
        self.ny = self.data.shape[0]

    def save_images(self):
        # self.raw_image.save(os.path.join(self.save_dir, self.raw_image_name))
        self.processed_image.save(os.path.join(self.save_dir, self.processed_image_name))

    def show_images(self):
        self.raw_image.show()
        self.processed_image.show()

    def _add_prefix(self, prefix, name, image_format):
        output_name = '{}_{}'.format(prefix, name) if prefix else name
        return '{}.{}'.format(output_name, image_format)

    def _check_files(self):
        if not os.path.isfile(self.file_path):
            raise ValueError('Provided file "{}" does not exist.'.format(self.file_path))

    def _check_input_type(self, file_path):
        self.possible_extensions = {
            'image': ['tif', 'tiff', 'png', 'bmp', 'gif', 'jpg', 'jpeg'],
            'npy': ['npy'],
        }
        extension = os.path.splitext(file_path)[1][1:].lower()
        for k in self.possible_extensions.keys():
            for e in self.possible_extensions[k]:
                if extension == e:
                    return k
        all_extensions = [x for x_list in self.possible_extensions.values() for x in x_list]
        all_extensions += [x.upper() for x in all_extensions]
        raise ValueError('Incorrect extension: {}. Possible values: {}.'.format(extension, ', '.join(all_extensions)))


# ********************** Create transmission element from the data from an image file:
def srwl_opt_setup_transm_from_file(
        file_path, resolution, thickness, delta, atten_len,
        arTr=None, extTr=0, fx=1e+23, fy=1e+23,
        xc=0, yc=0, ne=1, e_start=0, e_fin=0,
        area=None, rotate_angle=None, rotate_reshape=None,
        cutoff_background_noise=None, background_color=None,
        tile=None, shift_x=None, shift_y=None, invert=None,
        is_save_images=True, prefix='', output_image_format=None,
        _cutoff_max_fact=None, _max_color=255, _ret='srw' #OC28052020
):
    """Setup Sample element.

    :param file_path: path to the input file (image or .npy).
    :param resolution: resolution of the image [m/pixel].
    :param thickness: thickness of the sample [m].
    :param delta: refractive index decrement.
    :param atten_len: attenuation length [m].
    :param arTr: complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse coordinates.
    :param extTr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1).
    :param fx: estimated focal length in the horizontal plane [m].
    :param fy: estimated focal length in the vertical plane [m].
    :param xc: horizontal coordinate of center [m].
    :param yc: vertical coordinate of center [m].
    :param ne: number of transmission data points vs photon energy.
    :param e_start: initial photon energy [eV].
    :param e_fin: final photon energy [eV].
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_save_images: a flag to save the initial and processed images.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    :param _cutoff_max_fact: the ratio/factor for maximum signal cutoff (between 0 and 1).
    :param _max_color: the maximum color code to use for signal level above the cutoff_max_fact (0=black, 255=white).
    :param _ret: type of return: 'srw' returns SRWLOptT type transmission object which simulates the Sample (default), 'img' returns processed PIL image object, 'all' returns both SRWLOptT object and PIL image.
    :return: by default, transmission (SRWLOptT) type optical element which simulates the Sample; other options are defined by _ret variable.
    """

    input_parms = {
        "type": "sample",
        "resolution": resolution,
        "thickness": thickness,
        "refractiveIndex": delta,
        "attenuationLength": atten_len,
        "horizontalCenterCoordinate": xc,
        "verticalCenterCoordinate": yc,
        "initialPhotonEnergy": e_start,
        "finalPhotonPnergy": e_fin,
        'area': area,
        'rotateAngle': rotate_angle,
        'rotateReshape': rotate_reshape,
        'cutoffBackgroundNoise': cutoff_background_noise,
        'backgroundColor': background_color,
        'tile': tile,
        'shiftX': shift_x,
        'shiftY': shift_y,
        'invert': invert,
        'outputImageFormat': output_image_format,
        'cutoffMaxFact': _cutoff_max_fact,
        'maxColor': _max_color,
    }

    s = SRWLUtiSmp(
        file_path=file_path,
        area=area,
        rotate_angle=rotate_angle,
        rotate_reshape=rotate_reshape,
        cutoff_background_noise=cutoff_background_noise,
        background_color=background_color,
        tile=tile,
        shift_x=shift_x,
        shift_y=shift_y,
        invert=invert,
        is_show_images=False,
        is_save_images=is_save_images,
        prefix=prefix,
        output_image_format=output_image_format,
        cutoff_max_fact=_cutoff_max_fact, #OC29052020
        max_color=_max_color, #OC29052020
    )

    # Input parameters to SRWLOptT:
    nx = s.nx 
    ny = s.ny
    #print(nx, ny) #DEBUG
    
    rx = nx * resolution
    ry = ny * resolution

    #opT = srwlib.SRWLOptT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
    #                      _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
    #                      _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    data = s.data

    if(_ret == 'img'): return s.processed_image #OC28052020

    #OC10112018
    specPropAreDef = False
    if(ne > 1):
        if((isinstance(delta, list) or isinstance(delta, array)) and (isinstance(atten_len, list) or isinstance(atten_len, array))):
            lenDelta = len(delta)
            if((lenDelta == len(atten_len)) and (lenDelta == ne)): specPropAreDef = True
            else: raise Exception("Inconsistent spectral refractive index decrement and/or attenuation length data")

    #OC10112018
    useNumPy = False
    try:
        import numpy as np
        useNumPy = True
    except:
        print('NumPy can not be loaded, native Python arrays / lists will be used instead, impacting performance')

    if(useNumPy): #RC161018
        
        thickByLim = thickness/s.limit_value
        nxny = nx*ny
        miHalfThickByLimByAttenLen = None
        miThickDeltaByLim = None

        if(ne <= 1):
            miHalfThickByLimByAttenLen = -0.5*thickByLim/atten_len
            miThickDeltaByLim = -thickByLim*delta

            #amplTransm = np.exp(-0.5 * data * thickness / (s.limit_value * atten_len))
            amplTransm = np.exp(miHalfThickByLimByAttenLen*data)

            #optPathDiff =  -1 * data * thickness * delta / s.limit_value
            optPathDiff =  miThickDeltaByLim*data

            arTr = np.empty((2*nxny), dtype=float)

            #print(len(amplTransm), len(amplTransm[0])) #DEBUG
            
            arTr[0::2] = np.reshape(amplTransm, nxny)
            arTr[1::2] = np.reshape(optPathDiff, nxny)
            #opT.arTr = arTr
            
        else:
            two_ne = 2*ne
            arTr = np.empty((nxny*two_ne), dtype=float)

            if(specPropAreDef == False):
                miHalfThickByLimByAttenLen = -0.5*thickByLim/atten_len
                miThickDeltaByLim = -thickByLim*delta
            
            for ie in range(ne):
                if(specPropAreDef):
                    miHalfThickByLimByAttenLen = -0.5*thickByLim/atten_len[ie]
                    miThickDeltaByLim = -thickByLim*delta[ie]

                amplTransm = np.exp(miHalfThickByLimByAttenLen*data)
                optPathDiff = miThickDeltaByLim*data

                two_ie = 2*ie
                arTr[two_ie::two_ne] = np.reshape(amplTransm, nxny) #to check!
                arTr[(two_ie + 1)::two_ne] = np.reshape(optPathDiff, nxny)
    else:
        #Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
        nTot = 2*ne*nx*ny
        arTr = array('d', [0]*nTot)
    
        offset = 0
        for iy in range(ny):
            for ix in range(nx):
                #In images Y=0 corresponds from upper-left corner, in SRW it's lower-left corner:
                pathInBody = thickness * data[ny - iy - 1, ix] / s.limit_value
            
                #OC10112018
                miHalfPathInBody = -0.5*pathInBody
            
                if(specPropAreDef):
                    for ie in range(ne):
                        #opT.arTr[offset] = math.exp(-0.5 * pathInBody / atten_len)  # amplitude transmission
                        #opT.arTr[offset + 1] = -delta * pathInBody  # optical path difference
                        arTr[offset] = math.exp(miHalfPathInBody / atten_len[ie]) #amplitude transmission
                        arTr[offset + 1] = -delta[ie] * pathInBody #optical path difference
                        offset += 2
                else:
                    for ie in range(ne):
                        arTr[offset] = math.exp(miHalfPathInBody / atten_len) #amplitude transmission
                        arTr[offset + 1] = -delta * pathInBody  #optical path difference
                        offset += 2

    opT = srwlib.SRWLOptT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
                          _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
                          _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    opT.input_parms = input_parms

    #OC28052020
    if(_ret == 'srw'): return opT 
    elif(_ret == 'all'): return opT, s.processed_image
    else: return None
    #return opT

# ********************** Create transmission element from the data from a generated random 2D disk
def srwl_opt_setup_smp_rnd_obj2d( #RAC06052020
        _thickness, _delta, _atten_len, _rx = 10.e-06, _ry = 10.e-06, _xc = 0, _yc = 0, _nx = 1001, _ny = 1001,
        _dens = 1, _edge_frac = 0.02, _sim_step_size = 1.e-06,
        _obj_type = 1, _r_min_bw_obj = 0.1e-06, 
        _obj_size_min = 0.1e-06, _obj_size_max = 0.1e-06, _size_dist = 1, 
        _ang_min = 0, _ang_max = 0, _ang_dist = 1, _rand_alg = 1,
        _obj_par1 = None, _obj_par2 = None,
        _ext_tr = 1, _e_start = 0, _e_fin = 0, _ret = 'srw'): #OC01062020
        #_ext_tr = 1, _ne = 1, _e_start = 0, _e_fin = 0, _ret = 'srw'): #OC24052020
        #_extTr=1, _ne=1, _e_start=0, _e_fin=0, _ret='srw', _file_path='data_example_17', _file_name='ex17_res_opt_path_dif-ASCII'):
    """Setup Transmission type Optical Element which simulates Random 2D Disk using "srwl_uti_smp_rand_obj2d.py" 

    :param _thickness: thickness of the sample [m]
    :param _delta: refractive index decrement of sample material
    :param _atten_len: attenuation lengthof sample material [m]
    :param _rx: range of the horizontal coordinate [m] for which the transmission is defined
    :param _ry: range of the vertical coordinate [m] for which the transmission is defined
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :param _nx: number of transmission data points in the horizontal direction
    :param _ny: number of transmission data points in the vertical direction
    :param _dens: approximate density of objects [number of objects per mm^2]
    :param _edge_frac edge fraction defining how close to the edge can an object be placed on the object plane
    :param _sim_step_size: step size [m] that the simulation will take before trying to place an object (not used in some cases?)
    :param _obj_type: shape of objects to be generated: 1= rectangle, 2= ellipse, 3= triangle, 4= polygon, 5= random shapes
    :param _r_min_bw_obj:  minimum distance [m] between the centers of any two randomly placed objects 
    :param _obj_size_min: minimum allowed size of objects [m]
    :param _obj_size_max: maximum allowed size of objects [m]
    :param _size_dist: type of distribution of object sizes: 1- uniform, 2- normal (Gaussian), 3- Flory–Schulz
    :param _ang_min: minimum rotation angle of objects [degrees]
    :param _ang_max: maximum rotation angle of shape [degrees]
    :param _ang_dist: distribution of rotation angle. Choices are: 1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _rand_alg: randomization algorithm:
        1=  uniform seeding (default); algorithm will create a set of objects positioned on a uniform rectangular grid,
            based on density, minimum and maximum object size, and the minimum distance between objects;
            it then applies random seeded noise to the point locations;
        2=  2D random walk: each object will be positioned using a 2D random walk algorithm;
    :param _obj_par1: optional parameter defining object shapes; its meaning depends on the choice of object type (i.e. the value of _obj_type):
        if _obj_type == 1 (rectangle): it is ratio of width to length (length is the longest dimension);
            leaving it as None (default) will set the width to length ratio to 1, unless _obj_par2 = True
            is selected to randomize the width to length ratio;
        if _obj_type == 2 (ellipse): ratio of minor to major semi-axes;
            leaving value as None will set this ratio to 1, unless _obj_par2 = True
            is selected to randomize the minor to major semi-axes ratio;
        if _obj_type == 3 (triangle): ratio of height (y, or "opposite") to width (x, or adjunct);
            leaving it as None will set height to width ratio to 1 unless _obj_par2 = True
            is selected to randomize the height to width ratio (default is an equilateral triangle);
        if _obj_type == 4 (regular polygon): number of polygon vertices;
            leaving it as None will give a hexagon (6 sides) unless _obj_par2 = True is selected
            to randomize the numer of polygon vertices (max. 12);
        if _obj_type == 5 (random shapes): list of identifiers of shapes to randomly generate
            (1= square, 2= circle, 3= triangle, 4= hexagon, 5= any of 1-4);
            leaving it as None will generate all shape options as [1,2,3,4,5]
    :param _obj_par2: optional switch (True / False) defining randomization of object shapes, depending on object type (i.e. the value of _obj_type):
        if _obj_type == 1 (rectangle): randomizes the width to length ratio;
        if _obj_type == 2 (ellipse): randomizes the minor to major semi-axes ratio;
        if _obj_type == 3 (triangle): randomizes the height (y, or "opposite") to width (x, or adjunct) ratio;
        if _obj_type == 4 (regular polygon): randomizes the numer of polygon vertices;
    :param _extTr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1)
    :param _e_start: initial value of photon energy (active if > 0 and if _delta and _atten_len are arrays / lists)
    :param _e_fin: final value of photon energy (active if > 0 and if _delta and _atten_len are arrays / lists)
    :param _ret: type of return: 'srw' returns SRWLOptT type transmission object which simulates the Sample (default),
        'img' returns processed PIL image object, 'all' returns both SRWLOptT object and PIL image.
    :return: by default, transmission (SRWLOptT) type optical element which simulates the Sample; other options are defined by _ret variable
    """

    ### All input parameters to include in return object
    input_parms = { #what for?
        "type": "sample_rnd_obj2d",
        "thickness": _thickness,
        "refractiveIndex": _delta,
        "attenuationLength": _atten_len,
        "horizontalRange": _rx,
        "verticalRange": _ry,
        "horizontalCenterCoordinate": _xc,
        "verticalCenterCoordinate": _yc,
        "horizontalPoints": _nx,
        "verticalPoints": _ny,
        "density": _dens,
        "edgeFraction": _edge_frac,
        "simulationStepSize": _sim_step_size,
        "objectType": _obj_type,
        "minDistanceBetweenObjects": _r_min_bw_obj, 
        "objectMinSize": _obj_size_min,
        "objectMaxSize": _obj_size_max,
        "sizeDistributionType": _size_dist, 
        "minRotationAngle": _ang_min,
        "maxRotationAngle": _ang_max,
        "angleDistributionType": _ang_dist,
        "randomizationAlgorithm": _rand_alg,
        "objectParameter1": _obj_par1,
        "objectParameter2": _obj_par2,
        "externalTransmissionType": _ext_tr,
        #"pointNumberPhotonEnergy": _ne, 
        "initialPhotonEnergy": _e_start,
        "finalPhotonPnergy": _e_fin,}

    import_err_msg = '"{}" library cannot be imported. Please install it first with "pip install {}".'

    ### Load package Numpy
    try:
        import numpy as np
    except:
        raise ImportError(import_err_msg.format('NumPy', 'numpy')) #OC24052020
        #print('NumPy can not be loaded. You may need to install numpy. If you are using pip, you can use the following ' + 
        #      'command to install it: \npip install numpy')

    ### Load package srwl_uti_smp_rand_obj2d.py
    try:
        #sys.path.append(os.path.join(os.path.dirname(__file__), '..')) #OC24052020 (commented-out)
        from srwl_uti_smp_rnd_obj2d import get_rnd_2D, on_pxy, get_r1j, uni_rnd_seed
    except:
        raise Exception('srwl_uti_smp_rnd_obj2d module failed to load. Please make sure that skimage module is installed. It can be installed with "pip install scikit-image" ') #OC24052020
        #print('srwl_uti_smp_rnd_obj2d.py functions get_rnd_2D, on_pxy, get_r1j, and uni_rnd_seed failed to load.')

    ### Begin time record of creating and placing objects
    #t0 = time.time();

    ### Test User Input Parameters
    if _obj_size_min > _obj_size_max:
        raise Exception('Minimum object size (_obj_size_min) is greater than maximum object size (_obj_size_max).') #OC24052020
        #print('Minimum object size (_obj_size_min) is greater than maximum object size (_obj_size_max).')
    if _ang_min > _ang_max:
        raise Exception('Minimum rotation angle (_ang_min) is greater than maximum rotation angle (_ang_max).') #OC24052020  
        #print('Minimum rotation angle (_ang_min) is greater than maximum rotation angle (_ang_max).' )    

    ### Generate random 2D disk with random uniform positions (that will be random-walk moved later)
      #Initilize the object center positions
    rx_pixels = _nx #:param rx_pixels = _nx: (srwl_uti_smp_rand_obj2d.py) box size (x), namely each frame size [array size (pixels)]
    ry_pixels = _ny #:param ry_pixels = _ny: (srwl_uti_smp_rand_obj2d.py) box size (y), namely each frame size [array size (pixels)]
    num = int(_dens * ((_rx*1000)*(_ry*1000))) #:param num: (srwl_uti_smp_rand_obj2d.py) particle number.

    ### Generate random 2D disk (get_rnd_2D)
    if _rand_alg == 2: # 2D random walk
      #Simulate a twoD random walk 
      #Each of the 2D points, for each step, move with random direction and with amplitue as amp 
        #:param px/py: numpy.meshgrid, coordination of x/y, shape will be [point number, point number]
        #:param amp: float, each step move amptitde
        #:param rmin: the min distance between any two points
        #:param try_max: try number to acheive rmin, if can't achevie print that point
        #:return: new px/py
        px = np.random.uniform( _edge_frac, 1-_edge_frac, num ) * rx_pixels #numpy.random.uniform(low, high, size)
        py = np.random.uniform( _edge_frac, 1-_edge_frac, num ) * ry_pixels #numpy.random.uniform(low, high, size) 
        amp = (_sim_step_size * px.shape[0])/_rx #convert simulation_step_size from units [m] to [# of cells]
        rmin = (_r_min_bw_obj * px.shape[0])/_rx #convert minimum_obj_dist from units [m] to [# of cells]
        #try:
        px, py = get_rnd_2D( px, py, amp = amp, rmin = rmin, try_max=10000)
        #print("srwl_uti_smp_rnd_obj2d returned a numpy meshgrid of " + str(num) + " objects ... ") #OC24052020 (commented-out)
        #except:
        #    print("srwl_uti_smp_rnd_obj2d failed to return a numpy meshgrid")
    else: #_rand_alg == 1 # random ~uniform seeding (default)
      #Random ~uniform seeding (_rand_alg = 1 #default): Algorithm will 
      # create a uniform set of objects based on: density, minimum and 
      # maximum object size, and the minimum distance between objects. 
      # It then applies random seeded noise to the point locations.
        rmin = (_r_min_bw_obj * num)/_rx #convert minimum_obj_dist from units [m] to [# of cells]
        size_max = int(_obj_size_max * (((rx_pixels/(_rx*2))+(rx_pixels/(_rx*2)))/2))
        #try:
        px, py = uni_rnd_seed(num, rx_pixels, ry_pixels, obj_max_size=size_max, min_dist=rmin) 
        #print("srwl_uti_smp_rnd_obj2d returned a numpy meshgrid of " + str(num) + " objects ... ")
        #except:
        #    print("srwl_uti_smp_rnd_obj2d failed to return a numpy meshgrid")

    ### Get radius of particles in units of pixels (to be used with srwl_uti_smp_rand_obj2d.py def on_pxy)
      #particle radius: size[pixels] = particle_radius [m] * (((rx_pixels/(_rx [m] *2))+(rx_pixels/(_rx [m] *2)))/2)
    size_max = int(_obj_size_max * (((rx_pixels/(_rx*2))+(rx_pixels/(_rx*2)))/2))
    size_min = int(_obj_size_min * (((rx_pixels/(_rx*2))+(rx_pixels/(_rx*2)))/2))
    if (size_min < 1):
        raise Exception('Particle radius too small. Try using a larger particle minimum size.') #OC24052020  
        #print("\nParticle radius too small for \"srwl_uti_smp_rnd_obj2d.py\" function \"on_pxy\" to map onto disk." + 
        #      " Try using a larger particle minimum size.\n")

    ### Place the shape in position (on_pxy). px, py is numpy.meshgrid, shape will be [N,N]
    #try:
    object_plane = on_pxy(px, py, bx = rx_pixels, by = ry_pixels,
                          _obj_type = _obj_type, _obj_size_min = size_min, _obj_size_max = size_max, _size_dist = _size_dist,
                          _ang_min = _ang_min, _ang_max = _ang_max, _ang_dist = _ang_dist,
                          _obj_par1 = _obj_par1, _obj_par2 = _obj_par2)
    #print('completed (lasted', round(time.time() - t0, 6), 's)')
    #except:
    #print("\"srwl_uti_smp_rnd_obj2d.py\" function \"on_pxy\" failed to return a numpy array. " +
    #      "You may need to install numpy and scikit-image. If you are using pip, you can use the following command to install them: " + 
    #      "\npip install numpy scikit-image")

    ### Select output type: srw transmission object, png, jpg, or tiff
    ### Get data from processed array
    limit_value = 255
    # processed data
    data = ((np.subtract(object_plane, object_plane.min(), dtype=np.float32)) / 
            (np.subtract(object_plane.max(), object_plane.min(), dtype=np.float32)) * limit_value)

#OC24052020 (commented-out)
##    #Save plot as png
##    #Load package matplotlib
##    if _ret == 'png':
##        try:
##            import matplotlib as mpl
##        except:
##            print('Matplotlib can not be loaded. You may need to install matplotlib. If you are using pip, you can use the following ' + 
##              "command to install it: \npip install matplotlib")
##        filename = os.path.join(os.getcwd(), _file_path, _file_name)
##        return mpl.image.imsave(filename + '.png', data, cmap = 'binary')
##    
##    #Save plot as jpg
##    #Load package matplotlib
##    elif _ret == 'jpg':
##        try:
##            import matplotlib as mpl
##        except:
##            print('Matplotlib can not be loaded. You may need to install matplotlib. If you are using pip, you can use the following ' + 
##              "command to install it: \npip install matplotlib")
##        filename = os.path.join(os.getcwd(), _file_path, _file_name)
##        return mpl.image.imsave(filename + '.jpg', data, cmap = 'binary')
##    
##    #Save plot as tiff
##    #Load package matplotlib
##    elif _ret == 'tiff' or _ret == 'tif':
##        try:
##            import matplotlib
##            matplotlib.use('Agg')
##            import matplotlib.pyplot as plt
##        except:
##            print('Matplotlib can not be loaded. You may need to install matplotlib. If you are using pip, you can use the following ' + 
##              "command to install it: \npip install matplotlib")
##        filename = os.path.join(os.getcwd(), _file_path, _file_name)
##
##        dpi = 200 # Arbitrary. The number of pixels in the image will always be identical
##        tif_height, tif_width = np.array(data.shape, dtype=float) / dpi
##        fig = plt.figure(figsize=(tif_width, tif_height), dpi=dpi)
##        ax = fig.add_axes([0, 0, 1, 1])
##        ax.axis('off')
##        ax.imshow(data, interpolation='none', cmap = 'binary')
##        
##        return fig.savefig(filename + '.tif', dpi=dpi, cmap = 'binary')
##    
##    else: #_ret == 'srw'

    #OC24052020
    img = None
    if((_ret == 'img') or (_ret == 'all')): #Return of PIL image object is required

        try:
            from PIL import Image
        except ImportError:
            raise ImportError(import_err_msg.format('PIL', 'pillow'))

        img = Image.fromarray(data.astype(np.uint8))
        if(_ret == 'img'): return img
 
    #elif(_ret == 'srw' or _ret == 'all'): #OC28052020
    #elif _ret == 'srw':
    ### Create _void_cen_rad: flat array/list of void center coordinates and radii: [x1, y1, r1, x2, y2, r2,...]
    #OC24052020 (commented-out)
    #_void_cen_rad_X = np.ravel(px)
    #_void_cen_rad_Y = np.ravel(py)
    #_void_cen_rad_R = np.ones(len(_void_cen_rad_X)) * ((_obj_size_max + _obj_size_min)/2)
    #_void_cen_rad = [None]*(len(_void_cen_rad_X)+len(_void_cen_rad_Y)+len(_void_cen_rad_R))
    #_void_cen_rad[::3] = _void_cen_rad_X # add X coordinates (starting at postion 0 with step size 3)
    #_void_cen_rad[1::3] = _void_cen_rad_Y # add Y coordinates (starting at postion 1 with step size 3)
    #_void_cen_rad[2::3] = _void_cen_rad_R # add radius size (starting at postion 2 with step size 3)

    # spectral refractive index decrement and/or attenuation length data
    ne = 1 #OC01062020
    #specPropAreDef = False # is spectral propagation defined?
    #if(_ne > 1): # if number of transmission data points vs photon energy > 1

    deltaIsArray = isinstance(_delta, list) or isinstance(_delta, array) or isinstance(_delta, np.ndarray) #OC01062020
    attenLenIsArray = isinstance(_atten_len, list) or isinstance(_atten_len, array) or isinstance(_atten_len, np.ndarray)
    #if((isinstance(_delta, list) or isinstance(_delta, array) or isinstance(_delta, np.ndarray)) and 
    #   (isinstance(_atten_len, list) or isinstance(_atten_len, array) or isinstance(_atten_len, np.ndarray))):
        #    lenDelta = len(_delta) # length refractive index list
        #    if((lenDelta == len(_atten_len)) and (lenDelta == _ne)): specPropAreDef = True
        #else: raise Exception("Inconsistent spectral refractive index decrement and/or attenuation length data")
    if deltaIsArray and attenLenIsArray: #OC01062020
        ne = len(_delta)
        ne1 = len(_atten_len)
        if(ne > ne1): ne = ne1
    else:
        if deltaIsArray:
            ne = len(_delta)
            arAttenLen = [_atten_len]*ne
            _atten_len = arAttenLen
        elif attenLenIsArray:
            ne = len(_atten_len)
            arDelta = [_delta]*ne
            _delta = arDelta

    # process numpy array using physical properties such as thickness/attenuation-length/etc.
    thickByLim = _thickness/limit_value # thickness by limit
    nxny = _nx*_ny # transmission data points horizontal x vertical 
    miHalfThickByLimByAttenLen = None # initalize half thick limit by attenuation length
    miThickDeltaByLim = None # initalize thick refractive index by limit

    #creating arTr: (srwlib.SRWLOptT) complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude 
    if(ne <= 1): #OC01062020 # if number of transmission data points vs photon energy <= 1
    #if(_ne <= 1): # if number of transmission data points vs photon energy <= 1
        miHalfThickByLimByAttenLen = -0.5*thickByLim/_atten_len # (-1/2)*(thickness/limit_value)/(attenuation length)
        miThickDeltaByLim = -thickByLim*_delta # -(thickness/limit_value)*(refractive index)

        # amplitude transmission data
        amplTransm = np.exp(miHalfThickByLimByAttenLen*data)
        # optical path difference
        optPathDiff =  miThickDeltaByLim*data 

        #Create arTr: (srwlib.SRWLOptT) complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude 
        # transmission and optical path difference as function of transverse coordinates
        arTr = np.empty((2*nxny), dtype=float) 
        arTr[0::2] = np.reshape(amplTransm, nxny)
        arTr[1::2] = np.reshape(optPathDiff, nxny)

    else: # if number of transmission data points vs photon energy > 1
        #two_ne = 2*_ne
        two_ne = 2*ne #OC01062020
        arTr = np.empty((nxny*two_ne), dtype=float)

        #if(specPropAreDef == False):
        #    miHalfThickByLimByAttenLen = -0.5*thickByLim/_atten_len# (-1/2)*(thickness/limit_value)/(attenuation length)
        #    miThickDeltaByLim = -thickByLim*_delta# -(thickness/limit_value)*(refractive index)

        #for ie in range(_ne):
        for ie in range(ne): #OC01062020
            #if(specPropAreDef):
            #OC01062020
            miHalfThickByLimByAttenLen = -0.5*thickByLim/_atten_len[ie] # (-1/2)*(thickness/limit_value)/(attenuation length)[ie]
            miThickDeltaByLim = -thickByLim*_delta[ie] # -(thickness/limit_value)*(refractive index)[ie]

            # amplitude transmission data
            amplTransm = np.exp(miHalfThickByLimByAttenLen*data)
            # optical path difference
            optPathDiff = miThickDeltaByLim*data

            #Create arTr: (srwlib.SRWLOptT) complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse position
            # transmission and optical path difference as function of transverse coordinates
            two_ie = 2*ie
            arTr[two_ie::two_ne] = np.reshape(amplTransm, nxny)
            arTr[(two_ie + 1)::two_ne] = np.reshape(optPathDiff, nxny)

    ### Create transmission (SRWLOptT) type optical element
    #opT = srwlib.SRWLOptT(_nx=_nx, _ny=_ny, _rx=_rx, _ry=_ry, _arTr=arTr, _extTr=_ext_tr, _x=_xc, _y=_yc, _ne=_ne, _eStart=_e_start, _eFin=_e_fin)
    opT = srwlib.SRWLOptT(_nx=_nx, _ny=_ny, _rx=_rx, _ry=_ry, _arTr=arTr, _extTr=_ext_tr, _x=_xc, _y=_yc, _ne=ne, _eStart=_e_start, _eFin=_e_fin) #OC01062020
    opT.input_parms = input_parms # add user input parameters to transmission optical element

    #OC28052020
    if(_ret == 'srw'): return opT 
    elif(_ret == 'all'): return opT, img
    else: return None
