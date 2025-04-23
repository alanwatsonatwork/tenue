"""
This module configures tenue for 2k CCD on the Harlingten 50 cm telescope.
"""

import tenue.instrument

import numpy as np


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
)
