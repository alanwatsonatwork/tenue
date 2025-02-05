import numpy as np

import tenue.instrument


def overscanyslice(header):
    binning = int(header["SDTBN"])
    return slice(int(0 / binning), int(4096 / binning))


def overscanxslice(header):
    binning = int(header["SDTBN"])
    return slice(int(0 / binning), int(48 / binning))


def trimyslice(header):
    binning = int(header["SDTBN"])
    return slice(int(0 / binning), int(4096 / binning))


def trimxslice(header):
    binning = int(header["SDTBN"])
    return slice(int(48 / binning), int(4144 / binning))


def dorotate(header, data):
    rotation = header["SMTMRO"]
    return np.rot90(data, -int(rotation / 90))


def alpha(header):
    return header["SMTMRA"]


def delta(header):
    return header["SMTMDE"]


def rotation(header):
    date = header["DATE-OBS"][:10]
    return 0


def pixelscale(header):
    binning = int(header["SDTBN"])
    return 0.38 * binning / 3600


def gain(header):
    return 2.23


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
    gain=gain,
)
