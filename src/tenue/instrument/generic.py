"""
This module defines a generic instrument for Tenue.
"""

import tenue.instrument


def overscanyslice(header):
    return None


def overscanxslice(header):
    return None


def trimyslice(header):
    return None


def trimxslice(header):
    return None

def gain(header):
    return 1


tenue.instrument.setvalues(
    overscanxslice=overscanxslice,
    overscanyslice=overscanyslice,
    trimxslice=trimxslice,
    trimyslice=trimyslice,
    gain=gain,
)
