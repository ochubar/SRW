# -*- coding: utf-8 -*-
#############################################################################
# SRW Samples library
# Authors/Contributors: Maksim Rakitin
# October 26, 2016
#############################################################################

import math
import os

import srwlib
import uti_io
from srwlib import SRWLOpt


# ********************** The class for Samples:
class SRWLOptSmp(SRWLOpt):
    """The class for Samples from image file (.tif), NumPy array (.npy), etc.

    :param file_path: full path to the image or the saved NumPy array.
    :param bottom_limit: the bottom limit separating the image and the legend (black block).
    :param cutoff_background: the ratio for cutoff the background noise.
    :param is_show_images: a flag to show the initial and processed images.
    :param is_save_images: a flag to save the initial and processed images.
    :param raw_image_name: the name of the raw file in case if it's saved.
    :param processed_image_name: the name of the processed file in case if it's saved.
    :param prefix: the prefix to add to the names of the saved image files.
    """

    def __init__(self, file_path, bottom_limit=None, cutoff_background=0.5, is_show_images=False, is_save_images=False,
                 raw_image_name='raw', processed_image_name='processed', prefix='', output_image_format='tif'):
        # Input parameters:
        super(self.__class__, self).__init__()
        self.file_path = file_path
        self.bottom_limit = bottom_limit
        self.cutoff_background = cutoff_background
        self.is_show_images = is_show_images
        self.is_save_images = is_save_images
        self.raw_image_name = self._add_prefix(prefix, raw_image_name, output_image_format)
        self.processed_image_name = self._add_prefix(prefix, processed_image_name, output_image_format)

        # Aux variables:
        self._import_error_msg = '"{0}" library is not installed. Use "pip install {0}" to install it if you want to use this functionality.'

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
        d = uti_io.read_image(
            image_path=self.file_path,
            bottom_limit=self.bottom_limit,
            cutoff_background=self.cutoff_background,
        )
        self.data = d['data']
        self.raw_image = d['raw_image']
        self.processed_image = d['processed_image']
        self.bottom_limit = d['bottom_limit']
        self.limit_value = d['limit_value']

    def get_image_from_data(self):
        try:
            import numpy as np
        except:
            raise ImportError('The specified npy file cannot be read since {}'.format(self._import_error_msg.format('numpy')))
        data = np.load(self.file_path)
        self.limit_value = 255
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        try:
            from PIL import Image
            self.raw_image = Image.fromarray(self.data.astype(np.uint8))
            self.processed_image = self.raw_image
        except:
            print('The processed image will not be saved since {}'.format(self._import_error_msg.format('pillow')))
            self.is_save_images = False

    def read_sample(self):
        if self.input_type == 'image':
            self.get_data_from_image()
        elif self.input_type == 'npy':
            self.get_image_from_data()
        else:
            raise NotImplementedError(
                'Processing of the "{}" input type is not implemented yet.'.format(self.input_type))
        self.nx = self.data.shape[0]
        self.ny = self.data.shape[1]

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
def srwl_opt_setup_transmission_from_file(file_path, resolution, thickness, delta, atten_len,
                                          xc=0, yc=0, e_start=0, e_fin=0,
                                          is_save_images=True, prefix=''):
    """Setup Sample element.

    :param file_path: path to the input file (.tif or .npy).
    :param limit_value: maximum possible value (from bits per point).
    :param nx: number of horizontal points.
    :param ny: number of vertical points.
    :param resolution: resolution of the image [m/pixel].
    :param thickness: thickness of the sample [m].
    :param delta: refractive index decrement.
    :param atten_len: attenuation length [m].
    :param xc: horizontal coordinate of center [m].
    :param yc: vertical coordinate of center [m].
    :param e_start: initial photon energy [eV].
    :param e_fin: final photon energy [eV].
    :param is_save_images: a flag to save the initial and processed images.
    :param prefix: the prefix to add to the names of the saved image files.
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
    }

    s = SRWLOptSmp(
        file_path=file_path,
        is_show_images=False,
        is_save_images=is_save_images,
        prefix=prefix,
    )

    # Input parameters to SRWLOptT:
    nx = s.nx
    ny = s.ny
    rx = nx * resolution
    ry = ny * resolution
    fx = 1e+23
    fy = 1e+23
    ne = 1

    opT = srwlib.SRWLOptT(nx, ny, rx, ry, None, 1, fx, fy, xc, yc, ne, e_start, e_fin)

    # Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    offset = 0
    for iy in range(ny):
        for ix in range(nx):
            for ie in range(ne):
                pathInBody = thickness * s.data[ix, iy] / s.limit_value
                opT.arTr[offset] = math.exp(-0.5 * pathInBody / atten_len)  # amplitude transmission
                opT.arTr[offset + 1] = -delta * pathInBody  # optical path difference
                offset += 2

    opT.input_parms = input_parms

    return opT
