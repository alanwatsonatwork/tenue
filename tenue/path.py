import glob
import os.path

import tenue.fits
import tenue.instrument


def getrawfitspaths(fitspaths, filter=None, exposuretime=None):
    fitspaths = sorted(glob.glob(fitspaths))
    if filter != None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if tenue.instrument.filter(tenue.fits.readrawheader(fitspath)) == filter
        )
    if exposuretime != None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if tenue.instrument.exposuretime(tenue.fits.readrawheader(fitspath))
            == exposuretime
        )
    return fitspaths
