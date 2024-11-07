import glob
import os.path

import tenue.fits
import tenue.instrument


def getrawfitspaths(fitspaths, filter=None, exposuretime=None, fitspathsslice=None):
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
