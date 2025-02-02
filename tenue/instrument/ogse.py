import numpy as np

import tenue.instrument


def overscanyslice(header):
    binning = int(header["SDTBN"])
    return slice(int(2 / binning), int(20 / binning))


def overscanxslice(header):
    binning = int(header["SDTBN"])
    return slice(int(36 / binning), int(4132 / binning))


def trimyslice(header):
    binning = int(header["SDTBN"])
    return slice(int(22 / binning), int(4118 / binning))


def trimxslice(header):
    binning = int(header["SDTBN"])
    return slice(int(36 / binning), int(4132 / binning))


def dorotate(header, data):
    rotation = header["SMTMRO"]
    return np.rot90(data, -int(rotation / 90))


def alpha(header):
    return header["SMTMRA"]


def delta(header):
    return header["SMTMDE"]


def rotation(header):
    date = header["DATE-OBS"][:10]
    if date < "2025-01-01":
        return 90
    else:
        return 62


def pixelscale(header):
    binning = int(header["SDTBN"])
    return 0.20 * binning / 3600


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
