import numpy as np

import tenue.instrument


def overscanyslice(header):
    return slice(2, 20)


def overscanxslice(header):
    return slice(36, 4132)


def trimyslice(header):
    return slice(22, 4118)


def trimxslice(header):
    return slice(18, 4132)


def dorotate(header, data):
    rotation = header["SMTMRO"]
    return np.rot90(data, -int(rotation / 90))


def alpha(header):
    return header["SMTMRA"]


def delta(header):
    return header["SMTMDE"]


def rotation(header):
    return 90
    return 63.703


def pixelscale(header):
    return 0.20 / 3600


tenue.instrument.setvalues(
    overscanxslice=overscanxslice,
    overscanyslice=overscanyslice,
    trimxslice=trimxslice,
    trimyslice=trimyslice,
    dorotate=dorotate,
    rotation=rotation,
    alpha=alpha,
    delta=delta,
    pixelscale=pixelscale,
)
