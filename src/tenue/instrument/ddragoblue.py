import numpy as np
import math

import tenue.instrument
import tenue.image


def dooverscan(name, header, data):
    # Currently only works for the 4kx4k window with binning 1.
    binning = int(header["SDTBN"])
    if data.shape[0] == (4096 // binning) and data.shape[1] == (4192 // binning):
        p0 = tenue.image.clippedmean(
            data[(0 // binning) : (2048 // binning), (0 // binning) : (48 // binning)],
            sigma=3,
        )
        p1 = tenue.image.clippedmean(
            data[
                (0 // binning) : (2048 // binning),
                (4144 // binning) : (4192 // binning),
            ],
            sigma=3,
        )
        p2 = tenue.image.clippedmean(
            data[
                (2048 // binning) : (4096 // binning), (0 // binning) : (48 // binning)
            ],
            sigma=3,
        )
        p3 = tenue.image.clippedmean(
            data[
                (2048 // binning) : (4096 // binning),
                (4144 // binning) : (4192 // binning),
            ],
            sigma=3,
        )
        print(
            "%s: removing overscan levels of %.2f, %.2f, %.2f, and %.2f DN."
            % (name, p0, p1, p2, p3)
        )
        data[
            (0 // binning) : (2048 // binning), (0 // binning) : (2096 // binning)
        ] -= p0
        data[
            (0 // binning) : (2048 // binning), (2096 // binning) : (4192 // binning)
        ] -= p1
        data[
            (2048 // binning) : (4096 // binning), (0 // binning) : (2096 // binning)
        ] -= p2
        data[
            (2048 // binning) : (4096 // binning), (2096 // binning) : (4192 // binning)
        ] -= p3
    else:
        print("%s: skipping removing overscan levels from windowed data." % name)


def trimyslice(header):
    binning = int(header["SDTBN"])
    return slice(int(0 / binning), int(4096 / binning))


def trimxslice(header):
    binning = int(header["SDTBN"])
    return slice(int(48 / binning), int(4144 / binning))


def dorotate(header, data):
    if header["DTDS"] == "SI 1110-167":
        rotation = -header["SMTMRO"]
    elif header["DTDS"] == "SI 1110-185":
        rotation = header["SMTMRO"] + 90   
    return np.rot90(data, -int(rotation / 90))


def alpha(header):
    return header["STRSTRA"]


def delta(header):
    return header["STRSTDE"]


def rotation(header):
    date = header["DATE-OBS"][:10]
    return 0


def pixelscale(header):
    binning = int(header["SDTBN"])
    return 0.38 * binning / 3600
    
def _boresightparameters(header):
    if header["DTDS"] == "SI 1110-167":
        return -190, 0, math.radians(int(header["EMTMRO"]))
    elif header["DTDS"] == "SI 1110-185":
        return 110, 60, math.radians(int(header["EMTMRO"] + 90))

def boresightdx(header):
    boresightdx, boresightdy, theta = _boresightparameters(header)
    return boresightdx * math.cos(theta) - boresightdy * math.sin(theta)


def boresightdy(header):
    boresightdx, boresightdy, theta = _boresightparameters(header)
    return boresightdx * math.sin(theta) + boresightdy * math.cos(theta)


def gain(header):
    return 2.23


tenue.instrument.setvalues(
    dooverscan=dooverscan,
    trimxslice=trimxslice,
    trimyslice=trimyslice,
    dorotate=dorotate,
    rotation=rotation,
    alpha=alpha,
    delta=delta,
    pixelscale=pixelscale,
    boresightdx=boresightdx,
    boresightdy=boresightdy,
    gain=gain,
)
