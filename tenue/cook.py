import math
import os.path

import numpy as np

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
        raise RuntimeError("no bias found.")
    return _biasdata


def readdark(exposuretime, path="dark-{exposuretime:.0f}.fits", name="readdark"):
    global _darkdata
    path = path.format(exposuretime=exposuretime)
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _darkdata = tenue.fits.readproductdata(path)
    else:
        raise RuntimeError("no dark found.")
    return _darkdata


def readflat(filter, path="flat-{filter}.fits", name="readflat"):
    global _flatdata
    path = path.format(filter=filter)
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _flatdata = tenue.fits.readproductdata(path)
    else:
        raise RuntimeError("no flat found.")
    return _flatdata


def readmask(filter, path="mask-{filter}.fits", name="readmask"):
    global _maskdata
    path = path.format(filter=filter)
    if os.path.exists(path):
        print("%s: reading %s." % (name, path))
        _maskdata = tenue.fits.readproductdata(path)
    else:
        raise RuntimeError("no mask found.")
    return _maskdata


def writebias(path="bias.fits", name="writebias"):
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _biasdata)
    return


def writedark(path="dark-{exposuretime:.0f}.fits", exposuretime=None, name="writebias"):
    path = path.format(exposuretime=exposuretime)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _darkdata)
    return


def writeflat(path="flat-{filter}.fits", filter=None, name="writeflat"):
    path = path.format(filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _flatdata, filter=filter)
    return


def writemask(path="mask-{filter}.fits", filter=None, name="writemask"):
    path = path.format(filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _maskdata, filter=filter)
    return


def cook(
    fitspath,
    name="cook",
    dooverscan=False,
    dotrim=False,
    dobias=False,
    dodark=False,
    doflat=False,
    domask=False,
    dosky=False,
    dowindow=False,
    dorotate=False,
):

    print("%s: reading file %s." % (name, os.path.basename(fitspath)))
    header = tenue.fits.readrawheader(fitspath)
    data = tenue.fits.readrawdata(fitspath)

    # Set invalid pixels to nan.
    data[np.where(data == tenue.instrument.datamax(header))] = np.nan

    if dooverscan:
        overscandata = data[
            tenue.instrument.overscanyslice(header), tenue.instrument.overscanxslice(header)
        ]
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        sigma = tenue.image.clippedsigma(overscandata, sigma=3)
        print("%s: removing overscan level of %.2f ± %.2f DN." % (name, mean, sigma))
        data -= mean

    if dotrim:
        print("%s: trimming." % (name))
        data = data[tenue.instrument.trimyslice(header), tenue.instrument.trimxslice(header)]

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

    if dosky:
        median = np.nanmedian(data)
        print("%s: subtracting median sky of %.1f DN." % (name, median))
        # data -= np.nanmedian(data, axis=0, keepdims=True)
        # data -= np.nanmedian(data, axis=1, keepdims=True)
        data -= np.nanmedian(data, keepdims=True)

    if dorotate:
        print("%s: rotating to standard orientation." % (name))
        data = tenue.instrument.dorotate(header, data)

    if dowindow:
        print("%s: windowing." % (name))
        cx = 1996 / 2 + 30
        cy = 2028 / 2 - 0
        n = 512 * 3
        sx = int(cx - n / 2)
        sy = int(cy - n / 2)
        data = data[sy : sy + n, sx : sx + n]

    return data


def usefakebias():
    global _biasdata
    _biasdata = np.zeros((2051, 1024))
    return _biasdata


def usefakedark():
    global _darkdata
    _darkdata = np.zeros((2051, 1024))
    return _darkdata


def usefakeflat():
    global _flatdata
    _flatdata = np.ones((2051, 1024))
    return _flatdata


def usefakemask():
    global _maskdata
    _maskdata = np.ones((2051, 1024))
    return _maskdata


def makebias(fitspaths, biaspath="bias.fits"):
    def readonebias(fitspath):
        return cook(fitspath, name="makebias", dooverscan=True, dotrim=True)

    print("makebias: making bias from %s." % (fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(fitspaths)
    if len(fitspathlist) == 0:
        print("ERROR: no bias files found.")
        return

    datalist = list(readonebias(fitspath) for fitspath in fitspathlist)

    if len(datalist) == 0:
        print("ERROR: no bias files found.")
        return

    print("makebias: averaging %d biases with rejection." % len(datalist))
    global _biasdata
    _biasdata, biassigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = tenue.image.clippedmeanandsigma(_biasdata, sigma=5)
    print("makebias: bias is %.2f ± %.2f DN." % (mean, sigma))

    sigma = tenue.image.clippedmean(biassigma, sigma=5) / math.sqrt(len(datalist))
    print("makebias: estimated noise in bias is %.2f DN." % sigma)

    tenue.image.show(_biasdata, zscale=True)

    writebias(biaspath, name="makebias")

    print("makebias: finished.")

    return


def makedark(fitspaths, exposuretime, darkpath="dark-{exposuretime}.fits"):

    def readonedark(fitspath):
        return cook(
            fitspath, name="makedark", dooverscan=True, dotrim=True, dobias=True
        )

    print("makedark: making %.0f second dark from %s." % (exposuretime, fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(fitspaths, exposuretime=exposuretime)
    if len(fitspathlist) == 0:
        print("ERROR: no dark files found.")
        return

    datalist = list(readonedark(fitspath) for fitspath in fitspathlist)

    if len(datalist) == 0:
        print("ERROR: no dark files found.")
        return

    print("makedark: averaging %d darks with rejection." % len(datalist))
    global _darkdata
    _darkdata, darksigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    mean, sigma = tenue.image.clippedmeanandsigma(_darkdata, sigma=5)
    print("makedark: dark is %.2f ± %.2f DN." % (mean, sigma))

    sigma = tenue.image.clippedmean(darksigma, sigma=5) / math.sqrt(len(datalist))
    print("makedark: estimated noise in dark is %.2f DN." % sigma)

    tenue.image.show(_darkdata, zscale=True)

    writedark(darkpath, exposuretime=exposuretime, name="makedark")

    print("makedark: finished.")

    return


def makeflatandmask(
    fitspaths, filter, flatpath="flat-{filter}.fits", maskpath="mask-{filter}.fits"
):

    def readoneflat(fitspath):
        header = tenue.fits.readrawheader(fitspath)
        data = cook(
            fitspath,
            name="makeflatandmask",
            dooverscan=True,
            dotrim=True,
            dobias=True,
            dodark=True,
            domask=True,
        )
        if np.isnan(data).all():
            print("makeflatandmask: rejecting image: no valid data.")
            return None
        median = np.nanmedian(data)
        print("makeflatandmask: median is %.2f DN." % median)
        if median > tenue.instrument.flatmax(header):
            print("makeflatandmask: rejecting image: median too high.")
            return None
        else:
            print("makeflatandmask: accepting image.")
            data /= median
            return data

    def makeflathelper():
        datalist = list(readoneflat(fitspath) for fitspath in fitspathlist)
        datalist = list(data for data in datalist if data is not None)
        print("makeflatandmask: averaging %d flats with rejection." % (len(datalist)))
        flatdata, flatsigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

        mean, sigma = tenue.image.clippedmeanandsigma(flatdata, sigma=5)
        print("makeflatandmask: flat is %.2f ± %.3f." % (mean, sigma))

        sigma = tenue.image.clippedmean(flatsigma, sigma=5) / math.sqrt(len(datalist))
        print("makeflatandmask: estimated noise in flat is %.4f." % sigma)

        return flatdata

    def makemaskhelper(flatdata):

        maskdata = np.ones(flatdata.shape, dtype="float32")

        print("makeflatandmask: masking nan values.")
        maskdata[np.isnan(flatdata)] = 0

        print("makeflatandmask: masking inf values.")
        maskdata[np.isinf(flatdata)] = 0

        print("makeflatandmask: masking globally low pixels.")
        maskdata[np.where(flatdata < 0.80)] = 0

        print("makeflatandmask: masking locally high or low pixels.")
        low = tenue.image.medianfilter(flatdata, 7)
        high = flatdata / low
        maskdata[np.where(high < 0.97)] = 0
        maskdata[np.where(high > 1.03)] = 0

        print("makeflatandmask: masking pixels with at least two masked neighbors.")
        # Grow the mask so that any pixel with at least 2 neigboring bad pixels is also bad.
        grow = tenue.image.uniformfilter(maskdata, size=3)
        maskdata[np.where(grow <= 7 / 9)] = 0

        print(
            "makeflatandmask: fraction of masked pixels is %.4f."
            % (1 - np.nanmean(maskdata))
        )

        return maskdata

    print("makeflatandmask: making %s flat and mask from %s." % (filter, fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(fitspaths, filter=filter)
    if len(fitspathlist) == 0:
        print("ERROR: no flat files found.")
        return

    print("makeflatandmask: making fake mask.")
    usefakemask()

    print("makeflatandmask: making flat with fake mask.")
    flatdata = makeflathelper()

    print("makeflatandmask: making real mask.")
    global _maskdata
    _maskdata = makemaskhelper(flatdata)

    tenue.image.show(_maskdata, zrange=True)

    writemask(maskpath, filter=filter, name="makeflatandmask")

    print("makeflatandmask: making flat with real mask.")
    global _flatdata
    _flatdata = makeflathelper()

    tenue.image.show(_flatdata, zrange=True)

    writeflat(flatpath, filter=filter, name="makeflatandmask")

    print("makeflatandmask: finished.")

    return
