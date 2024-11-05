import glob
import os.path

import tenue.fits
import tenue.instrument


def getrawfitspaths(directorypath, filter=None, exposuretime=None):
    fitspaths = glob.glob(directorypath + "/*.fits")
    if len(fitspaths) == 0:
        fitspaths = glob.glob(directorypath + "/*.fits.fz")
    fitspaths = sorted(fitspaths)
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
            if tenue.fits.readrawheader(fitspath)[tenue.instrument.exposuretimekeyword()]
            == exposuretime
        )
    return fitspaths
