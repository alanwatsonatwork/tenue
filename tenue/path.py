import glob
import os.path

from tenue.fits import readrawheader

def getrawfitspaths(directorypath, filter=None):
    fitspaths = sorted(glob.glob(directorypath + "/*.fits"))
    if filter == None:
        return fitspaths
    else:
        return list(
            fitspath
            for fitspath in fitspaths
            if readrawheader(fitspath)["FILTER"] == filter
        )
