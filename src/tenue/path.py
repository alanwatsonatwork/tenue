import glob
import os.path

import tenue.fits
import tenue.instrument


def getrawfitspaths(fitspaths, filter=None, exposuretime=None, fitspathsslice=None):
    """
    Return an expanded and filtered list of FITS paths.

    The ``fitspath`` argument will be expanded by :func:`glob.glob`. This expansion must give
    a list of names of FITS files or compressed FITS files. If ``filter`` is not
    ``None``, then files which do not have that filter value are eliminated
    from the list. If ``exposuretime`` is not ``None``, then files which do not
    have that exposure time are eliminated from the list. Finally, if
    ``fitspathsslice`` is not ``None``, then this list is sliced.

    :param fitspaths: A patter to be expanded by :func:`glob.glob`. The
        expansion must give a list of names of FITS files or compressed FITS
        files.

    :param filter: A filter name as a string or ``None``. If it is not ``None``,
        the files which do not have that exposure time are        eliminated
        from the expanded list.
    
    :param exposuretime: An exposure time as a number or ``None``. If it is not
        ``None``, then files which do not have that exposure time are eliminated
        from the expanded list.
    
    :param fitspathsslice: A slice, either of the string ``"firsthalf"`` or
        ``"secondhalf"``, or ``None``. If it is not ``None``, then the expanded
        and filtered list is sliced before being returned. The special values
        ``"firsthalf"`` or ``"secondhalf"`` refer, as might be expected, to
        slices containing the first half and second half of the file names.
    
    :return: An expanded, filtered, and sliced list of FITS file name.
    """
    fitspaths = sorted(glob.glob(fitspaths))
    if filter is not None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if tenue.instrument.filter(tenue.fits.readrawheader(fitspath)) == filter
        )
    if exposuretime is not None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if tenue.instrument.exposuretime(tenue.fits.readrawheader(fitspath))
            == exposuretime
        )
    if fitspathsslice == "firsthalf":
        fitspathsslice = slice(None, len(fitspaths) // 2)
    elif fitspathsslice == "secondhalf":
        fitspathsslice = slice(len(fitspaths) // 2, None)
    if fitspathsslice is not None:
        fitspaths = fitspaths[fitspathsslice]
    return fitspaths
