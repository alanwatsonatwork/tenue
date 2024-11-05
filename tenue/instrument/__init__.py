def setvalues(
    datamax=None,
    overscanyslice=None,
    overscanxslice=None,
    trimyslice=None,
    trimxslice=None,
    filterkeyword=None,
    exposuretimekeyword=None,
    alphakeyword=None,
    deltakeyword=None,
    rotationkeyword=None,
    pixelscale=None,
    rotation=None,
    rotationpositive=None,
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

    global _filterkeyword
    _filterkeyword = filterkeyword

    global _exposuretimekeyword
    _exposuretimekeyword = exposuretimekeyword

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

    global _rotationpositive
    _rotationpositive = rotationpositive

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


def filterkeyword():
    return _filterkeyword


def exposuretimekeyword():
    return _exposuretimekeyword


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


def rotationpositive():
    return _rotationpositive


def flatmax():
    return _flatmax
