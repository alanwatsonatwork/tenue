from datetime import datetime


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
    return datetime.fromisoformat(header["DATE-OBS"] + "Z").timestamp()


def defaultendtimestamp(header):
    return starttimestamp(header) + exposuretime(header)


def defaultrotation(header):
    return header["CROTA2"]


def defaultalpha(header):
    return header["CRVAL1"]


def defaultdelta(heafer):
    return header["CRVAL2"]


def defaultpixelscale(header):
    return header["CDELT2"]


def defaultgain(header):
    return header["GAIN"]


def setvalues(
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
    gain=defaultgain,
):

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

    global _gain
    _gain = gain


def datamax(header):
    return _datamax(header)


def flatmax(header):
    return _flatmax(header)


def overscanyslice(header):
    return _overscanyslice(header)


def overscanxslice(header):
    return _overscanxslice(header)


def trimyslice(header):
    return _trimyslice(header)


def trimxslice(header):
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


def gain(header):
    return _gain(header)
