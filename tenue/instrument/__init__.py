def setvalues(
    datamax=None,
    overscanyslice=None,
    overscanxslice=None,
    trimyslice=None,
    trimxslice=None,
    filterkeyword=None,
    alphakeyword=None,
    deltakeyword=None,
    pixelscale=None,
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

    global _filterkeyword
    _filterkeyword = filterkeyword

    global _alphakeyword
    _alphakeyword = alphakeyword

    global _deltakeyword
    _deltakeyword = deltakeyword

    global _pixelscale
    _pixelscale = pixelscale


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


def filterkeyword():
    return _filterkeyword


def alphakeyword():
    return _alphakeyword


def deltakeyword():
    return _deltakeyword


def pixelscale():
    return _pixelscale


print("instrument")
