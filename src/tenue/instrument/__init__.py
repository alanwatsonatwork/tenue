from datetime import datetime

import tenue.image
import math


def defaultdooverscan(name, header, data):
    if overscanyslice(header) is None or overscanxslice(header) is None:
        return
    overscandata = data[
        overscanyslice(header),
        overscanxslice(header),
    ]
    mean = tenue.image.clippedmean(overscandata, sigma=3)
    mean = tenue.image.clippedmean(overscandata, sigma=3)
    mean = tenue.image.clippedmean(overscandata, sigma=3)
    sigma = tenue.image.clippedsigma(overscandata, sigma=3)
    print("%s: removing overscan level of %.2f ± %.2f DN." % (name, mean, sigma))
    data -= mean


def defaultdatamax(header):
    return 65535


def defaultflatmax(header):
    return 32767


def defaultdorotate(header, data):
    return data


def defaultfilter(header):
    return header["FILTER"]


def defaultexposuretime(header):
    return header["EXPTIME"]


def defaultstarttimestamp(header):
    return datetime.fromisoformat(header["DATE-OBS"] + "+00:00").timestamp()


def defaultendtimestamp(header):
    return starttimestamp(header) + exposuretime(header)


def defaultrotation(header):
    return header["CROTA2"]


def defaultalpha(header):
    return header["CRVAL1"]


def defaultdelta(header):
    return header["CRVAL2"]


def defaultpixelscale(header):
    return header["CDELT2"]


def defaultboresightdx(header):
    return 0


def defaultboresightdy(header):
    return 0


def defaultgain(header):
    return header["GAIN"]


def setvalues(
    dooverscan=defaultdooverscan,
    datamax=defaultdatamax,
    flatmax=defaultflatmax,
    overscanyslice=None,
    overscanxslice=None,
    trimyslice=None,
    trimxslice=None,
    dorotate=defaultdorotate,
    filter=defaultfilter,
    exposuretime=defaultexposuretime,
    starttimestamp=defaultstarttimestamp,
    endtimestamp=defaultendtimestamp,
    rotation=defaultrotation,
    alpha=defaultalpha,
    delta=defaultdelta,
    pixelscale=defaultpixelscale,
    boresightdx=defaultboresightdx,
    boresightdy=defaultboresightdy,
    gain=defaultgain,
):

    global _dooverscan
    _dooverscan = dooverscan

    global _datamax
    _datamax = datamax

    global _flatmax
    _flatmax = flatmax

    global _overscanyslice
    _overscanyslice = overscanyslice

    global _overscanxslice
    _overscanxslice = overscanxslice

    global _trimyslice
    _trimyslice = trimyslice

    global _trimxslice
    _trimxslice = trimxslice

    global _dorotate
    _dorotate = dorotate

    global _filter
    _filter = filter

    global _exposuretime
    _exposuretime = exposuretime

    global _starttimestamp
    _starttimestamp = starttimestamp

    global _endtimestamp
    _endtimestamp = endtimestamp

    global _rotation
    _rotation = rotation

    global _alpha
    _alpha = alpha

    global _delta
    _delta = delta

    global _pixelscale
    _pixelscale = pixelscale

    global _boresightdx
    _boresightdx = boresightdx

    global _boresightdy
    _boresightdy = boresightdy

    global _gain
    _gain = gain


def dooverscan(name, header, data):
    return _dooverscan(name, header, data)


def datamax(header):
    return _datamax(header)


def flatmax(header):
    return _flatmax(header)


def overscanyslice(header):
    if _overscanyslice is None:
        return None
    else:
        return _overscanyslice(header)


def overscanxslice(header):
    if _overscanxslice is None:
        return None
    else:
        return _overscanxslice(header)


def trimyslice(header):
    if _trimyslice is None:
        return None
    else:
        return _trimyslice(header)


def trimxslice(header):
    if _trimxslice is None:
        return None
    else:
        return _trimxslice(header)


def dorotate(header, data):
    return _dorotate(header, data)


def filter(header):
    return _filter(header)


def exposuretime(header):
    return _exposuretime(header)


def starttimestamp(header):
    return _starttimestamp(header)


def endtimestamp(header):
    return _endtimestamp(header)


def rotation(header):
    return _rotation(header)


def alpha(header):
    return _alpha(header)


def delta(header):
    return _delta(header)


def pixelscale(header):
    return _pixelscale(header)


def boresightdx(header):
    return _boresightdx(header)


def boresightdy(header):
    return _boresightdy(header)


def gain(header):
    return _gain(header)
