import math
import os.path

import numpy as np

import scipy.ndimage

import tenue.fits
import tenue.image
import tenue.instrument
import tenue.path

_biasdata = None
_darkdata = None
_flatdata = None
_maskdata = None


def readbias(path="bias.fits", name="readbias"):
    global _biasdata
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _biasdata = tenue.fits.readproductdata(path)
    else:
        print("%s: WARNING: no bias found; using a fake bias." % (name))
        makefakebias()
    return _biasdata


def readdark(path="dark.fits", name="readdark"):
    global _darkdata
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _darkdata = tenue.fits.readproductdata(path)
    else:
        print("%s: WARNING: no dark found; using a fake dark." % (name))
        makefakedark()
    return _darkdata


def readflat(filter=None, path="flat-%s.fits", name="readflat"):
    global _flatdata
    if filter is not None and "%s" in path:
        path = path % filter
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _flatdata = tenue.fits.readproductdata(path)
    else:
        print("%s: WARNING: no flat found; using a fake flat." % (name))
        makefakeflat()
    return _flatdata


def readmask(filter=None, path="mask-%s.fits", name="readmask"):
    global _maskdata
    if filter is not None and "%s" in path:
        path = path % filter
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _maskdata = tenue.fits.readproductdata(path)
    else:
        print("%s: WARNING: no mask found; using a fake mask." % (name))
        makefakemask()
    return _maskdata


def writebias(data, path="bias.fits", name="writebias"):
    print("%s: writing %s." % (name, path))
    global _biasdata
    _biasdata = data
    tenue.fits.writeproduct(path, data)
    return


def writedark(data, path="dark.fits", name="writebias"):
    print("%s: writing %s." % (name, path))
    global _darkdata
    _darkdata = data
    tenue.fits.writeproduct(path, data)
    return


def writeflat(data, path="flat.fits", filter=None, name="writeflat"):
    print("%s: writing %s." % (name, path))
    global _flatdata
    _flatdata = data
    tenue.fits.writeproduct(path, data, filter=filter)
    return


def writemask(data, path="mask.fits", filter=None, name="writemask"):
    print("%s: writing %s." % (name, path))
    global _maskdata
    _maskdata = data
    tenue.fits.writeproduct(path, data, filter=filter)
    return


def cook(
    fitspath,
    name="cook",
    dooverscan=False,
    dotrim=False,
    dobias=False,
    dodark=False,
    doflat=False,
    donormalize=False,
    domask=False,
    dosky=False,
    dowindow=False,
    dorotate=False,
):

    print("%s: reading file %s." % (name, os.path.basename(fitspath)))
    header = tenue.fits.readrawheader(fitspath)
    data = tenue.fits.readrawdata(fitspath)

    # Set invalid pixels to nan.
    data[np.where(data == tenue.instrument.datamax())] = np.nan

    if dooverscan:
        overscandata = data[
            tenue.instrument.overscanyslice(), tenue.instrument.overscanxslice()
        ]
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        sigma = tenue.image.clippedsigma(overscandata, sigma=3)
        print("%s: removing overscan level of %.2f ± %.2f DN." % (name, mean, sigma))
        data -= mean

    if dotrim:
        print("%s: trimming." % (name))
        data = data[tenue.instrument.trimyslice(), tenue.instrument.trimxslice()]

    if dobias:
        print("%s: subtracting bias." % (name))
        data -= _biasdata

    if dodark:
        print("%s: subtracting dark." % (name))
        data -= _darkdata

    if doflat:
        print("%s: dividing by flat." % (name))
        data /= _flatdata

    if domask:
        print("%s: masking." % (name))
        data[np.where(_maskdata == 0)] = np.nan

    median = np.nanmedian(data)
    print("%s: median in sample region is %.1f DN." % (name, median))

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
        if tenue.instrument.rotationpositive():
            data = np.rot90(data, +int(rotation / 90))
        else:
            data = np.rot90(data, -int(rotation / 90))

    if dowindow:
        print("%s: windowing." % (name))
        cx = 1996 / 2 + 30
        cy = 2028 / 2 - 0
        n = 512 * 3
        sx = int(cx - n / 2)
        sy = int(cy - n / 2)
        data = data[sy : sy + n, sx : sx + n]

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

    print("makebias: making bias from %s." % (directorypath))

    fitspathlist = tenue.path.getrawfitspaths(directorypath)
    if len(fitspathlist) == 0:
        print("ERROR: no bias files found.")
        return

    datalist = list(readonebias(fitspath) for fitspath in fitspathlist)

    if len(datalist) == 0:
        print("ERROR: no bias files found.")
        return

    print("makebias: averaging %d biases with rejection." % len(datalist))
    biasdata, biassigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = tenue.image.clippedmeanandsigma(biasdata, sigma=5)
    print("makebias: bias is %.2f ± %.2f DN." % (mean, sigma))

    sigma = tenue.image.clippedmean(biassigma, sigma=5) / math.sqrt(len(datalist))
    print("makebias: estimated noise in bias is %.2f DN." % sigma)

    tenue.image.show(biasdata, zscale=True)

    writebias(biasdata, name="makebias")

    print("makebias: finished.")

    return


def makedark(directorypath):
    def readonedark(fitspath):
        return cook(
            fitspath, name="makedark", dooverscan=True, dotrim=True, dobias=True
        )

    print("makedark: making dark from %s." % (directorypath))

    fitspathlist = tenue.path.getrawfitspaths(directorypath)
    if len(fitspathlist) == 0:
        print("ERROR: no dark files found.")
        return

    datalist = list(readonedark(fitspath) for fitspath in fitspathlist)

    if len(datalist) == 0:
        print("ERROR: no dark files found.")
        return

    print("makedark: averaging %d darks with rejection." % len(datalist))
    darkdata, darksigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = tenue.image.clippedmeanandsigma(darkdata, sigma=5)
    print("makedark: dark is %.2f ± %.2f DN." % (mean, sigma))

    sigma = tenue.image.clippedmean(darksigma, sigma=5) / math.sqrt(len(datalist))
    print("makedark: estimated noise in dark is %.2f DN." % sigma)

    tenue.image.show(darkdata, zscale=True)

    writedark(darkdata, name="makedark")

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
        flatdata, flatsigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

        mean, sigma = tenue.image.clippedmeanandsigma(flatdata, sigma=5)
        print("makeflatandmask: flat is %.2f ± %.2f DN." % (mean, sigma))

        sigma = tenue.image.clippedmean(flatsigma, sigma=5) / math.sqrt(len(datalist))
        print("makeflatandmask: estimated noise in flat is %.3f." % sigma)

        return flatdata

    def makemaskhelper(flatdata):

        maskdata = np.ones(flatdata.shape)

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

    print("makeflatandmask: making %s flat and mask from %s." % (filter, directorypath))

    fitspathlist = tenue.path.getrawfitspaths(directorypath, filter=filter)
    if len(fitspathlist) == 0:
        print("ERROR: no flat files found.")
        return

    print("makeflatandmask: making fake mask.")
    makefakemask()

    print("makeflatandmask: making flat with fake mask.")
    flatdata = makeflathelper()

    print("makeflatandmask: making real mask.")
    maskdata = makemaskhelper(flatdata)
    tenue.image.show(maskdata, zrange=True)

    writemask(maskdata, "mask-%s.fits" % filter, filter=filter, name="makeflatandmask")

    print("makeflatandmask: making flat with real mask.")
    flatdata = makeflathelper()

    writeflat(flatdata, "flat-%s.fits" % filter, filter=filter, name="makeflatandmask")

    tenue.image.show(flatdata, zrange=True)

    print("makeflatandmask: finished.")

    return
