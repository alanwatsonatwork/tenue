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
            if tenue.fits.readrawheader(fitspath)[tenue.instrument.filterkeyword()]
            == filter
        )
    if exposuretime != None:
        fitspaths = list(
            fitspath
            for fitspath in fitspaths
            if tenue.fits.readrawheader(fitspath)[
                tenue.instrument.exposuretimekeyword()
            ]
            == exposuretime
        )
    return fitspaths
