import numpy as np
import astropy.io.fits

ihdu = 0


def _ihdu(fitspath):
    """
    Return the HDU index for the actual data in a FITS file. This is 0
    for an uncompressed FITS file and 1 for a fpack-compressed FITS
    file. The suffix of the fitspath is used to determine if the file is
    compressed or not; if the suffix is ".fz", the file is assumed to be
    compressed.
    """
    if fitspath[-3:] == ".fz":
        return 1
    else:
        return 0


def readrawheader(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading header from raw file %s." % (name, os.path.basename(fitspath))
        )
    hdu = astropy.io.fits.open(fitspath)
    return hdu[_ihdu(fitspath)].header


def readrawdata(fitspath, name=None):
    if name is not None:
        print("%s: reading data from raw file %s." % (name, os.path.basename(fitspath)))
    hdu = astropy.io.fits.open(fitspath)
    return np.array(hdu[_ihdu(fitspath)].data, dtype=np.float32)

