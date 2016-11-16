# -*- coding: utf-8 -*-
#############################################################################
# SRW Samples library
# Authors/Contributors: Maksim Rakitin
# October 26, 2016
#############################################################################

import os
import math
import srwlib
import numpy as np
from PIL import Image

#********************** The class for Samples:
class SRWLSamples:
    """The class for Samples from image file (.tif), NumPy array (.npy), etc.

    :param file_path: full path to the image or the saved NumPy array.
    :param input_type: the type of the input ('image', 'npy' or 'txt').
    :param bottom_limit: the bottom limit separating the image and the legend (black block).
    :param cutoff_background: the ratio for cutoff the background noise.
    :param is_show_images: a flag to show the initial and processed images.
    :param is_save_images: a flag to save the initial and processed images.
    :param raw_image_name: the name of the raw file in case if it's saved.
    :param processed_image_name: the name of the processed file in case if it's saved.
    """
    def __init__(self, file_path, input_type,
                 bottom_limit=None, cutoff_background=0.5, is_show_images=False, is_save_images=False,
                 raw_image_name='raw_image.tif', processed_image_name='processed_image.tif'):
        # Input parameters:
        self.file_path = file_path
        self.input_type = input_type
        self.bottom_limit = bottom_limit
        self.cutoff_background = cutoff_background
        self.is_show_images = is_show_images
        self.is_save_images = is_save_images
        self.raw_image_name = raw_image_name
        self.processed_image_name = processed_image_name

        # Output parameters:
        self.data = None
        self.raw_image = None
        self.processed_image = None
        self.nx = None
        self.ny = None
        self.limit_value = None

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
        # Read the image:
        self.raw_image = Image.open(self.file_path)

        # Convert it to NumPy array:
        data = np.array(self.raw_image)

        # Get bits per point:
        mode_to_bpp = {'1': 1, 'L': 8, 'P': 8, 'I;16': 16, 'RGB': 24, 'RGBA': 32, 'CMYK': 32, 'YCbCr': 24, 'I': 32, 'F': 32}
        bpp = mode_to_bpp[self.raw_image.mode]
        self.limit_value = float(2 ** bpp - 1)

        # Get the bottom limit if it's not provided:
        if not self.bottom_limit:
            self.bottom_limit = np.where(data[:, 0] == 0)[0][0]

        # Remove the bottom black area:
        if data[:, 0].max() > 0:  # do not remove if the background is already black
            data = np.copy(data[:self.bottom_limit, :])
        self.data = np.transpose(data)

        # Remove background noise:
        idxs_less = np.where(self.data < self.limit_value * self.cutoff_background)
        self.data[idxs_less] = np.uint16(0)

        # Generate new image object to track the changes:
        self.processed_image = Image.fromarray(np.transpose(self.data))

    def get_image_from_data(self):
        data = self.read_data_from_npy()
        self.limit_value = 255
        self.data = (data - data.min()) / (data.max() - data.min()) * self.limit_value
        self.raw_image = Image.fromarray(self.data.astype(np.uint8))
        self.processed_image = self.raw_image

    def read_data_from_npy(self):
        return np.load(self.file_path)

    def read_sample(self):
        possible_values = ('image', 'npy', 'txt')
        if self.input_type == 'image':
            self.get_data_from_image()
        elif self.input_type == 'npy':
            self.get_image_from_data()
        elif self.input_type == 'txt':
            raise NotImplementedError('Import of "txt" transmission objects is not implemented yet.')
        else:
            raise ValueError('Incorrect input type: {}. Possible values: {}.'.format(self.input_type, ', '.join(possible_values)))

        self.nx = self.data.shape[0]
        self.ny = self.data.shape[1]

    def save_images(self):
        self.raw_image.save(self.raw_image_name)
        self.processed_image.save(self.processed_image_name)

    def show_images(self):
        self.raw_image.show()
        self.processed_image.show()

    def _check_files(self):
        if not os.path.isfile(self.file_path):
            raise ValueError('Provided file "{}" does not exist.'.format(self.file_path))


#********************** Create transmission element from the data from an image file:
def srwl_opt_setup_transmission_from_image(file_path, resolution, thickness, delta, atten_len,
                                           xc=0, yc=0, e_start=0, e_fin=0, input_type='image'):
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
    :param input_type: the type of the input ('image', 'npy' or 'txt').
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

    s = SRWLSamples(file_path, input_type=input_type, is_show_images=False, is_save_images=True)

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
