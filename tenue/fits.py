import numpy as np
import astropy.io.fits

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


def readraw(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading header and data from FITS file %s."
            % (name, os.path.basename(fitspath))
        )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    header = hdu[ihdu].header
    data = np.array(hdu[ihdu].data, dtype=np.float32)
    hdu.close()
    return header, data


def readrawheader(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading header from raw FITS file %s."
            % (name, os.path.basename(fitspath))
        )
    header, data = readraw(fitspath, None)
    return header


def readrawdata(fitspath, name=None):
    if name is not None:
        print("%s: reading data from raw file %s." % (name, os.path.basename(fitspath)))
    header, data = readraw(fitspath, None)
    return data


def readproduct(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading header from product file %s."
            % (name, os.path.basename(fitspath))
        )
    hdu = astropy.io.fits.open(fitspath)
    ihdu = _ihdu(fitspath)
    header = hdu[ihdu].header
    data = np.array(hdu[ihdu].data, dtype=np.float32)
    hdu.close()
    return header, data


def readproductheader(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading header from product file %s."
            % (name, os.path.basename(fitspath))
        )
    header, data = readproduct(fitspath, None)
    return header


def readproductdata(fitspath, name=None):
    if name is not None:
        print(
            "%s: reading data from product file %s."
            % (name, os.path.basename(fitspath))
        )
    header, data = readproduct(fitspath, None)
    return data


def writeproduct(fitspath, data, name=None, filter=None):
    if name is not None:
        print("%s: writing product file %s." % (name, os.path.basename(fitspath)))
    header = astropy.io.fits.Header()
    if filter is not None:
        header.append(("FILTER", filter))
    astropy.io.fits.writeto(fitspath, data, header, overwrite=True)
    return
