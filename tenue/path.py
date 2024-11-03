import glob
import os.path

import tenue.fits
import tenue.instrument


def getrawfitspaths(directorypath, filter=None):
    fitspaths = glob.glob(directorypath + "/*.fits")
    if len(fitspaths) == 0:
        fitspaths = glob.glob(directorypath + "/*.fits.fz")
    fitspaths = sorted(fitspaths)
    if filter == None:
        return fitspaths
    else:
        return list(
            fitspath
            for fitspath in fitspaths
            if tenue.fits.readrawheader(fitspath)[tenue.instrument.filterkeyword()]
            == filter
        )
