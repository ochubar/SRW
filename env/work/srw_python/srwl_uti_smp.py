# -*- coding: utf-8 -*-
#############################################################################
# SRW Samples library
# Authors/Contributors: Maksim Rakitin
# October 26, 2016
# October 27, 2017
#############################################################################

import math
import os

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
                 raw_image_name='raw', processed_image_name='processed', prefix='', output_image_format=None):
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
        assert 0 <= self.background_color <= 255, 'Background color ({}) should be between 0 and 255'.format(
            self.background_color)
        assert 0 <= self.cutoff_background_noise <= 1, 'Cutoff background noise ({}) should be between 0 and 1'.format(
            self.cutoff_background_noise)
        self.data[np.where(self.data < self.limit_value * self.cutoff_background_noise)] = np.uint16(
            self.background_color)

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
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
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
        area=None, rotate_angle=None, rotate_reshape=None, cutoff_background_noise=None,
        background_color=None, tile=None, shift_x=None, shift_y=None, invert=None,
        is_save_images=True, prefix='', output_image_format=None,
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
    :return: transmission (SRWLOptT) type optical element which simulates the Sample.
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
    )

    # Input parameters to SRWLOptT:
    nx = s.nx
    ny = s.ny
    rx = nx * resolution
    ry = ny * resolution

    opT = srwlib.SRWLOptT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
                          _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
                          _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    data = s.data

    # Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    offset = 0
    for iy in range(ny):
        for ix in range(nx):
            for ie in range(ne):
                # In images Y=0 corresponds from upper-left corner, in SRW it's lower-left corner:
                pathInBody = thickness * data[ny - iy - 1, ix] / s.limit_value
                opT.arTr[offset] = math.exp(-0.5 * pathInBody / atten_len)  # amplitude transmission
                opT.arTr[offset + 1] = -delta * pathInBody  # optical path difference
                offset += 2

    opT.input_parms = input_parms

    return opT
