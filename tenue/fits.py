import numpy as np
import astropy.io.fits

ihdu = 0


def readrawheader(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading header from raw file %s." % (name, os.path.basename(fitspath))
        )
    hdu = astropy.io.fits.open(fitspath)
    return hdu[ihdu].header


def readrawdata(fitspath, name=None):
    if name is not None:
        print("%s: reading data from raw file %s." % (name, os.path.basename(fitspath)))
    hdu = astropy.io.fits.open(fitspath)
    return np.array(hdu[ihdu].data, dtype=np.float32)

