from datetime import datetime


def defaultexposuretime(header):
    return header["EXPTIME"]


def defaultstarttimestamp(header):
    return datetime.fromisoformat(header["DATE-OBS"] + "Z").timestamp()


def defaultendtimestamp(header):
    return starttimestamp(header) + exposuretime(header)


def defaultdorotate(header, data):
    return data


def setvalues(
    datamax=None,
    overscanyslice=None,
    overscanxslice=None,
    trimyslice=None,
    trimxslice=None,
    dorotate=defaultdorotate,
    filterkeyword=None,
    exposuretime=defaultexposuretime,
    starttimestamp=defaultstarttimestamp,
    endtimestamp=defaultendtimestamp,
    alphakeyword=None,
    deltakeyword=None,
    rotationkeyword=None,
    pixelscale=None,
    rotation=None,
    flatmax=None,
):

    global _datamax
    _datamax = datamax

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

    global _filterkeyword
    _filterkeyword = filterkeyword

    global _exposuretime
    _exposuretime = exposuretime

    global _starttimestamp
    _starttimestamp = starttimestamp

    global _endtimestamp
    _endtimestamp = endtimestamp

    global _alphakeyword
    _alphakeyword = alphakeyword

    global _deltakeyword
    _deltakeyword = deltakeyword

    global _rotationkeyword
    _rotationkeyword = rotationkeyword

    global _pixelscale
    _pixelscale = pixelscale

    global _rotation
    _rotation = rotation

    global _flatmax
    _flatmax = flatmax


def datamax():
    return _datamax


def overscanyslice():
    return _overscanyslice


def overscanxslice():
    return _overscanxslice


def trimyslice():
    return _trimyslice


def trimxslice():
    return _trimxslice


def dorotate(header, data):
    return _dorotate(header, data)


def filterkeyword():
    return _filterkeyword


def exposuretime(header):
    return _exposuretime(header)


def starttimestamp(header):
    return _starttimestamp(header)


def endtimestamp(header):
    return _endtimestamp(header)


def alphakeyword():
    return _alphakeyword


def deltakeyword():
    return _deltakeyword


def rotationkeyword():
    return _rotationkeyword


def pixelscale():
    return _pixelscale


def rotation():
    return _rotation


def flatmax():
    return _flatmax
