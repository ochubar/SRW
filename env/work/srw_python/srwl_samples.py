# -*- coding: utf-8 -*-
#############################################################################
# SRW Samples library
# Authors/Contributors: Maksim Rakitin
# October 26, 2016
#############################################################################

import math
import srwlib

#********************** Create transmission element from the data from an image file:
def srwl_opt_setup_transmission_from_image(image_data, limit_value, nx, ny, resolution, thickness, delta, atten_len,
                                           xc=0, yc=0, e_start=0, e_fin=0):
    """Setup Sample element.

    :param image_data: data from the provided TIFF file.
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
    :return: transmission (SRWLOptT) type optical element which simulates the Sample.
    """

    input_parms = {
        "type": "sample",
        "horizontalPoints": nx,
        "verticalPoints": ny,
        "resolution": resolution,
        "thickness": thickness,
        "refractiveIndex": delta,
        "attenuationLength": atten_len,
        "horizontalCenterCoordinate": xc,
        "verticalCenterCoordinate": yc,
        "initialPhotonEnergy": e_start,
        "finalPhotonPnergy": e_fin,
    }

    ne = 1
    fx = 1e+23
    fy = 1e+23
    rx = nx * resolution
    ry = ny * resolution

    opT = srwlib.SRWLOptT(nx, ny, rx, ry, None, 1, fx, fy, xc, yc, ne, e_start, e_fin)

    # Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    ofst = 0
    for iy in range(ny):
        for ix in range(nx):
            for ie in range(ne):
                pathInBody = thickness * image_data[ix, iy] / limit_value
                opT.arTr[ofst] = math.exp(-0.5 * pathInBody / atten_len)  # amplitude transmission
                opT.arTr[ofst + 1] = -delta * pathInBody  # optical path difference
                ofst += 2

    opT.input_parms = input_parms

    return opT
