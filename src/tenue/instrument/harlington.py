"""
This module configures tenue for 2k CCD on the Harlingten 50-cm telescope.
"""

import tenue.instrument

import numpy as np


def alpha(header):
    component = header["OBJCTRA"].split(" ")
    alpha = 15.0 * (
        float(component[0]) + float(component[1]) / 60.0 + float(component[2]) / 3600.0
    )
    return alpha


def delta(header):
    component = header["OBJCTDEC"].split(" ")
    if component[0][0] == "-":
        delta = (
            float(component[0])
            - float(component[1]) / 60.0
            - float(component[2]) / 3600.0
        )
    else:
        delta = (
            float(component[0])
            + float(component[1]) / 60.0
            + float(component[2]) / 3600.0
        )
    return delta

def pixelscale(header):
    binning = int(2048 / header["NAXIS1"])
    return 0.80 * binning / 3600

def rotation(header):
    return 0.0

def overscanyslice(header):
    return None


def overscanxslice(header):
    return None


def trimyslice(header):
    return None


def trimxslice(header):
    return None


def dorotate(header, data):
    return np.flip(data, axis=1)


def gain(header):
    return 1


tenue.instrument.setvalues(
    overscanxslice=overscanxslice,
    overscanyslice=overscanyslice,
    trimxslice=trimxslice,
    trimyslice=trimyslice,
    dorotate=dorotate,
    gain=gain,
    alpha=alpha,
    delta=delta,
    pixelscale=pixelscale,
    rotation=rotation,
)
