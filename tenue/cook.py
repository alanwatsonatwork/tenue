import os.path
import numpy as np
import astropy.stats
import scipy.ndimage
import matplotlib.pyplot as plt


import tenue.instrument
from tenue.fits import (
    readrawdata,
    readrawheader,
    readproductdata,
    readproductheader,
    writeproduct,
)
from tenue.path import getrawfitspaths


def readbias(directorypath, name="readbias"):
    global _biasdata
    path = directorypath + "/bias.fits"
    if os.path.exists(path):
        print("%s: reading bias.fits" % (name))
        _biasdata = readproductdata(path)
    else:
        print("%s: WARNING: no bias found; using a fake bias." % (name))
        makefakebias()
    return _biasdata


def readdark(directorypath, name="readdark"):
    global _darkdata
    path = directorypath + "/dark.fits"
    if os.path.exists(path):
        print("%s: reading dark.fits" % (name))
        _darkdata = readproductdata(path)
    else:
        print("%s: WARNING: no dark found; using a fake dark." % (name))
        makefakedark()
    return _darkdata


def readflat(directorypath, filter, name="readflat"):
    global _flatdata
    path = directorypath + ("/flat-%s.fits" % filter)
    if os.path.exists(path):
        print("%s: reading flat-%s.fits" % (name, filter))
        _flatdata = readproductdata(path)
    else:
        print("%s: WARNING: no flat found; using a fake flat." % (name))
        makefakeflat()
    return _flatdata


def readmask(directorypath, filter, name="readmask"):
    global _maskdata
    path = directorypath + ("/mask-%s.fits" % filter)
    if os.path.exists(path):
        print("%s: reading mask-%s.fits" % (name, filter))
        _maskdata = readproductdata(path)
    else:
        print("%s: WARNING: no mask found; using a fake mask." % (name))
        makefakemask()
    return _maskdata


def writebias(directorypath, data, name="writebias"):
    print("%s: writing bias.fits" % (name))
    global _biasdata
    _biasdata = data
    writeproduct(directorypath + "/bias.fits", data)
    return


def writedark(directorypath, data, name="writebias"):
    print("%s: writing dark.fits" % (name))
    global _darkdata
    _darkdata = data
    writeproduct(directorypath + "/dark.fits", data)
    return


def writeflat(directorypath, data, filter, name="writeflat"):
    print("%s: writing flat-%s.fits" % (name, filter))
    global _flatdata
    _flatdata = data
    writeproduct(directorypath + ("/flat-%s.fits" % filter), data, filter=filter)
    return


def writemask(directorypath, data, filter, name="writemask"):
    print("%s: writing mask-%s.fits" % (name, filter))
    global _maskdata
    _maskdata = data
    writeproduct(directorypath + ("/mask-%s.fits" % filter), data, filter=filter)
    return


def cook(
    fitspath,
    name="cook",
    dooverscan=False,
    dotrim=False,
    dobias=False,
    doflat=False,
    donormalize=False,
    domask=False,
    dosky=False,
    dowindow=False,
    dorotate=False,
):

    print("%s: reading file %s." % (name, os.path.basename(fitspath)))
    header = readrawheader(fitspath)
    data = readrawdata(fitspath)

    # Set saturated pixels to nan.
    data[np.where(data == (2**16 - 1))] = np.nan

    if dooverscan:
        # Converted from BIASSEC.
        overscandata = data[
            tenue.instrument.overscanyslice(), tenue.instrument.overscanxslice()
        ]
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            overscandata, sigma=3, cenfunc="median", stdfunc="mad_std"
        )
        print("%s: removing overscan level of %.2f ± %.2f DN." % (name, mean, sigma))
        data -= mean

    if dotrim:
        # Converted from DATASEC.
        print("%s: trimming." % (name))
        data = data[tenue.instrument.trimyslice(), tenue.instrument.trimxslice()]

    if dobias:
        print("%s: subtracting bias." % (name))
        data -= _biasdata

    if doflat:
        print("%s: dividing by flat." % (name))
        data /= _flatdata

    if domask:
        print("%s: masking." % (name))
        data[np.where(_maskdata == 0)] = np.nan

    if dowindow:
        print("%s: windowing." % (name))
        cx = 1996 / 2 + 30
        cy = 2028 / 2 - 0
        n = 512 * 3
        sx = int(cx - n / 2)
        sy = int(cy - n / 2)
        data = data[sy : sy + n, sx : sx + n]

    median = np.nanmedian(data)
    print("%s: median is %.1f DN." % (name, median))

    if dosky:
        print("%s: subtracting sky." % (name))
        # data -= np.nanmedian(data, axis=0, keepdims=True)
        # data -= np.nanmedian(data, axis=1, keepdims=True)
        data -= np.nanmedian(data, keepdims=True)

    if donormalize:
        print("%s: normalizing to median." % (name))
        data /= median

    if dorotate:
        print("%s: rotating to standard orientation." % (name))
        rotation = header["SMTMRO"]
        data = np.rot90(data, -int(rotation / 90))

    return data


def makefakebias():
    global _biasdata
    _biasdata = np.zeros((2051, 1024))
    return _biasdata


def makefakedark():
    global _darkdata
    _darkdata = np.zeros((2051, 1024))
    return _darkdata


def makefakeflat():
    global _flatdata
    _flatdata = np.ones((2051, 1024))
    return _flatdata


def makefakemask():
    global _maskdata
    _maskdata = np.ones((2051, 1024))
    return _maskdata


def makebias(directorypath):
    def readonebias(fitspath):
        return cook(fitspath, name="makebias", dooverscan=True, dotrim=True)

    print("makebias: making bias.fits from %s." % (directorypath))

    fitspathlist = getrawfitspaths(directorypath + "/bias/")
    if len(fitspathlist) == 0:
        print("ERROR: no bias files found.")
        return

    datalist = list(readonebias(fitspath) for fitspath in fitspathlist)

    if len(datalist) == 0:
        print("ERROR: no bias files found.")
        return

    print("makebias: averaging %d biases with rejection." % len(datalist))
    mean, median, sigma = astropy.stats.sigma_clipped_stats(
        datalist, sigma=3, axis=0, cenfunc="median", stdfunc="mad_std"
    )
    biasdata = mean

    mean, median, sigma = astropy.stats.sigma_clipped_stats(
        biasdata, sigma=5, cenfunc="median", stdfunc="mad_std"
    )
    print("makebias: final residual bias level is %.2f ± %.2f DN." % (mean, sigma))

    print("makebias: plotting median of columns.")
    plt.figure()
    plt.plot(range(biasdata.shape[1]), np.nanmedian(biasdata, axis=0))
    plt.show()
    print("makebias: plotting median of rows.")
    plt.figure()
    plt.plot(range(biasdata.shape[0]), np.nanmedian(biasdata, axis=1))
    plt.show()

    writebias(directorypath, biasdata, name="makebias")

    print("makebias: finished.")

    return


def makedark(directorypath):
    def readonedark(fitspath):
        return cook(
            fitspath, name="makedark", dooverscan=True, dotrim=True, dobias=True
        )

    print("makedark: making dark.fits from %s." % (directorypath))

    fitspathlist = getrawfitspaths(directorypath + "/dark/")
    if len(fitspathlist) == 0:
        print("ERROR: no dark files found.")
        return

    datalist = list(readonedark(fitspath) for fitspath in fitspathlist)

    if len(datalist) == 0:
        print("ERROR: no dark files found.")
        return

    print("makedark: averaging %d darks with rejection." % len(datalist))
    mean, median, sigma = astropy.stats.sigma_clipped_stats(
        datalist, sigma=3, axis=0, cenfunc="median", stdfunc="mad_std"
    )
    darkdata = mean

    mean, median, sigma = astropy.stats.sigma_clipped_stats(
        darkdata, sigma=5, cenfunc="median", stdfunc="mad_std"
    )
    print("makedark: final residual dark level is %.2f ± %.2f DN." % (mean, sigma))

    print("makedark: plotting median of columns.")
    plt.figure()
    plt.plot(range(darkdata.shape[1]), np.nanmedian(darkdata, axis=0))
    plt.show()
    print("makedark: plotting median of rows.")
    plt.figure()
    plt.plot(range(darkdata.shape[0]), np.nanmedian(darkdata, axis=1))
    plt.show()

    writedark(directorypath, darkdata, name="makebias")

    print("makedark: finished.")

    return


def makeflatandmask(directorypath, filter):
    def readoneflat(fitspath):
        return cook(
            fitspath,
            name="makeflatandmask",
            dooverscan=True,
            dotrim=True,
            dobias=True,
            domask=True,
            donormalize=True,
        )

    def makeflathelper():
        datalist = list(readoneflat(fitspath) for fitspath in fitspathlist)
        print("makeflatandmask: averaging %d flats with rejection." % (len(datalist)))
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            datalist, sigma=3, axis=0, cenfunc="median", stdfunc="mad_std"
        )
        flatdata = mean
        return flatdata

    def makemaskhelper(flatdata):

        maskdata = np.ones(flatdata.shape)

        print("makeflatandmask: masking hot pixels.")
        maskdata[np.where(_darkdata > 200)] = 0

        print("makeflatandmask: masking nan values.")
        maskdata[np.isnan(flatdata)] = 0

        print("makeflatandmask: masking inf values.")
        maskdata[np.isinf(flatdata)] = 0

        print("makeflatandmask: masking globally low pixels.")
        maskdata[np.where(flatdata < 0.80)] = 0

        print("makeflatandmask: masking locally high or low pixels.")
        low = scipy.ndimage.median_filter(flatdata, 7)
        high = flatdata / low
        maskdata[np.where(high < 0.97)] = 0
        maskdata[np.where(high > 1.03)] = 0

        print("makeflatandmask: masking pixels with at least two masked neighbors.")
        # Grow the mask so that any pixel with at least 2 neigboring bad pixels is also bad.
        grow = scipy.ndimage.filters.uniform_filter(maskdata, size=3, mode="nearest")
        maskdata[np.where(grow <= 7 / 9)] = 0

        print(
            "makeflatandmask: fraction of masked pixels is %.4f."
            % (1 - np.nanmean(maskdata))
        )

        return maskdata

    print("makeflatandmask: making %s flat from %s." % (filter, directorypath))

    fitspathlist = getrawfitspaths(directorypath + "/flat/", filter=filter)
    if len(fitspathlist) == 0:
        print("ERROR: no flat files found.")
        return

    readbias(directorypath, name="makeflatandmask")

    print("makeflatandmask: making fake mask.")
    makefakemask()

    print("makeflatandmask: making flat with fake mask.")
    flatdata = makeflathelper()

    print("makeflatandmask: making real mask.")
    maskdata = makemaskhelper(flatdata)

    writemask(directorypath, maskdata, filter, name="makeflatandmask")

    readmask(directorypath, filter, name="makeflatandmask")

    print("makeflatandmask: making flat with real mask.")
    flatdata = makeflathelper()

    writeflat(directorypath, flatdata, filter, name="makeflatandmask")

    print("makeflatandmask: plotting median of columns.")
    plt.figure()
    plt.plot(range(flatdata.shape[1]), np.nanmedian(flatdata, axis=0))
    plt.show()
    print("makeflatandmask: plotting median of rows.")
    plt.figure()
    plt.plot(range(flatdata.shape[0]), np.nanmedian(flatdata, axis=1))
    plt.show()

    print("makeflatandmask: finished.")

    return
