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


def writebias(path="bias.fits", name="writebias"):
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _biasdata, exposuretime=0)
    return


def writedark(path="dark-{exposuretime:.0f}.fits", exposuretime=None, name="writebias"):
    path = path.format(exposuretime=exposuretime)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _darkdata, exposuretime=exposuretime)
    return


def writeflat(path="flat-{filter}.fits", filter=None, name="writeflat"):
    path = path.format(filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _flatdata, filter=filter)
    return


def cook(
    fitspath,
    name="cook",
    dooverscan=False,
    dotrim=False,
    dobias=False,
    dodark=False,
    doflat=False,
    dosky=False,
    dowindow=False,
    dorotate=False,
):

    print("%s: reading %s." % (name, os.path.basename(fitspath)))
    header, data = tenue.fits.readraw(fitspath)

    # Set invalid pixels to nan.
    data[np.where(data == tenue.instrument.datamax(header))] = np.nan

    if dooverscan:
        overscandata = data[
            tenue.instrument.overscanyslice(header),
            tenue.instrument.overscanxslice(header),
        ]
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        mean = tenue.image.clippedmean(overscandata, sigma=3)
        sigma = tenue.image.clippedsigma(overscandata, sigma=3)
        print("%s: removing overscan level of %.2f ± %.2f DN." % (name, mean, sigma))
        data -= mean

    if dotrim:
        print("%s: trimming." % (name))
        data = data[
            tenue.instrument.trimyslice(header), tenue.instrument.trimxslice(header)
        ]

    if dobias:
        print("%s: subtracting bias." % (name))
        data -= _biasdata

    if dodark:
        print("%s: subtracting dark." % (name))
        data -= _darkdata

    if doflat:
        print("%s: dividing by flat." % (name))
        data /= _flatdata

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

    return header, data


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


def makebias(fitspaths, biaspath="bias.fits", fitspathsslice=None):

    print("makebias: making bias from %s." % (fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(
        fitspaths, exposuretime=0, fitspathsslice=fitspathsslice
    )

    if len(fitspathlist) == 0:
        print("ERROR: no bias files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = cook(fitspath, name="makebias", dooverscan=True, dotrim=True)
        headerlist.append(header)
        datalist.append(data)

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


def makedark(
    fitspaths, exposuretime, darkpath="dark-{exposuretime}.fits", fitspathsslice=None
):

    print("makedark: making %.0f second dark from %s." % (exposuretime, fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(
        fitspaths, exposuretime=exposuretime, fitspathsslice=fitspathsslice
    )

    if len(fitspathlist) == 0:
        print("ERROR: no dark files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = cook(
            fitspath, name="makedark", dooverscan=True, dotrim=True, dobias=True
        )
        headerlist.append(header)
        datalist.append(data)

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


def makeflat(fitspaths, filter, flatpath="flat-{filter}.fits", fitspathsslice=None):

    ############################################################################

    print("makeflat: making %s flat %s." % (filter, fitspaths))

    ############################################################################

    print("makeflat: making flat without mask.")

    fitspathlist = tenue.path.getrawfitspaths(
        fitspaths, filter=filter, fitspathsslice=fitspathsslice
    )

    if len(fitspathlist) == 0:
        print("ERROR: no flat files found.")
        return

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = cook(
            fitspath,
            name="makeflat",
            dooverscan=True,
            dotrim=True,
            dobias=True,
            dodark=True,
        )
        if np.isnan(data).all():
            print("makedark: rejected %s: no valid data." % os.path.basename(fitspath))
            continue
        median = np.nanmedian(data)
        print("makeflat: median is %.2f DN." % median)
        if median > tenue.instrument.flatmax(header):
            print("makeflat: rejecting image: median too high.")
            continue
        print("makeflat: accepted %s." % os.path.basename(fitspath))
        data /= median
        headerlist.append(header)
        datalist.append(data)

    print("makeflat: averaging %d flats with rejection." % (len(datalist)))

    flatdata, flatsigma = tenue.image.clippedmeanandsigma(datalist, sigma=3, axis=0)

    ############################################################################

    print("makeflat: making mask.")

    maskdata = np.ones(flatdata.shape, dtype="float32")

    print("makeflat: masking nan values.")
    maskdata[np.isnan(flatdata)] = 0

    print("makeflat: masking inf values.")
    maskdata[np.isinf(flatdata)] = 0

    print("makeflat: masking globally low pixels.")
    maskdata[np.where(flatdata < 0.80)] = 0

    print("makeflat: masking locally high or low pixels.")
    low = tenue.image.medianfilter(flatdata, 7)
    high = flatdata / low
    maskdata[np.where(high < 0.97)] = 0
    maskdata[np.where(high > 1.03)] = 0

    print("makeflat: masking pixels with at least two masked neighbors.")
    # Grow the mask so that any pixel with at least 2 neigboring bad pixels is also bad.
    grow = tenue.image.uniformfilter(maskdata, size=3)
    maskdata[np.where(grow <= 7 / 9)] = 0

    print("makeflat: fraction of masked pixels is %.4f." % (1 - np.nanmean(maskdata)))

    tenue.image.show(maskdata, zrange=True)

    ############################################################################

    print("makeflat: making flat with mask.")

    maskeddatalist = []
    for data in datalist:
        data[np.where(maskdata == 0)] = np.nan
        data / np.nanmedian(data)
        maskeddatalist.append(data)

    print("makeflat: averaging %d flats with rejection." % (len(maskeddatalist)))
    flatdata, flatsigma = tenue.image.clippedmeanandsigma(
        maskeddatalist, sigma=3, axis=0
    )

    mean, sigma = tenue.image.clippedmeanandsigma(flatdata, sigma=5)
    print("makeflat: flat is %.2f ± %.3f." % (mean, sigma))

    sigma = tenue.image.clippedmean(flatsigma, sigma=5) / math.sqrt(len(maskeddatalist))
    print("makeflat: estimated noise in flat is %.4f." % sigma)

    global _flatdata
    _flatdata = flatdata
    tenue.image.show(_flatdata, zrange=True)
    writeflat(flatpath, filter=filter, name="makeflat")

    ############################################################################

    print("makeflat: finished.")

    return
