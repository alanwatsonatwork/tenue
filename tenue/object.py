import math
import os.path
import numpy as np
import matplotlib.pyplot as plt

import tenue.cook
import tenue.fits
import tenue.image
import tenue.instrument
import tenue.path

_skydata = None

def writeobject(directorypath, data, filter, name="writeobject"):
    print("%s: writing object-%s.fits" % (name, filter))
    tenue.fits.writeproduct(
        directorypath + ("/object-%s.fits" % filter), data, filter=filter
    )
    return

def writesky(directorypath, data, filter, name="writesky"):
    print("%s: writing sky-%s.fits" % (name, filter))
    global _skydata
    _skydata = data
    tenue.fits.writeproduct(
        directorypath + ("/sky-%s.fits" % filter), data, filter=filter
    )
    return


def makeobject(
    directorypath,
    filter,
    align=None,
    nalignregion=40,
    refalpha=None,
    refdelta=None,
    sigma=None,
    nwindow=None,
    showalignment=True,
    doskyimage=False,
):
    def readonepointing(fitspath):
        header = tenue.fits.readrawheader(fitspath)
        print(
            "makeobject: reading pointing for %s object file %s."
            % (filter, os.path.basename(fitspath))
        )
        alpha = math.radians(header[tenue.instrument.alphakeyword()])
        delta = math.radians(header[tenue.instrument.deltakeyword()])
        print(
            "makeobject: pointing is alpha = %.5f deg delta = %.5f deg."
            % (math.degrees(alpha), math.degrees(delta))
        )
        return [alpha, delta]

    def readonesky(fitspath):

        print(
            "makeobject: reading %s sky file %s." % (filter, os.path.basename(fitspath))
        )

        data = tenue.cook.cook(
            fitspath,
            name="makeobject",
            dooverscan=True,
            dotrim=True,
            dobias=True,
            dodark=True,
            doflat=True,
            domask=True,
            dowindow=False,
            dosky=True,
            dorotate=True,
        )
        return data

    def readoneobject(fitspath):

        print(
            "makeobject: reading %s object file %s."
            % (filter, os.path.basename(fitspath))
        )

        data = tenue.cook.cook(
            fitspath,
            name="makeobject",
            dooverscan=True,
            dotrim=True,
            dobias=True,
            dodark=True,
            doflat=True,
            domask=True,
            dowindow=False,
            dosky=True,
            dorotate=True,
        )
        data -= _skydata

        header = tenue.fits.readrawheader(fitspath)

        alpha = math.radians(header[tenue.instrument.alphakeyword()])
        delta = math.radians(header[tenue.instrument.deltakeyword()])
        pixelscale = math.radians(tenue.instrument.pixelscale())
        rotation = math.radians(tenue.instrument.rotation())

        print(
            "makeobject: pointing is alpha = %.5f deg delta = %.5f deg."
            % (math.degrees(alpha), math.degrees(delta))
        )
        dalpha = (alpha - refalpha) / pixelscale * math.cos(refdelta)
        ddelta = (delta - refdelta) / pixelscale
        dx = int(np.round(dalpha * math.cos(rotation) - ddelta * math.sin(rotation)))
        dy = -int(np.round(dalpha * math.sin(rotation) + ddelta * math.cos(rotation)))
        print("makeobject: raw offset is dx = %+d px dy = %+d px." % (dx, dy))

        margin = 512
        if align != None:
            aligny = align[0] - margin
            alignx = align[1] - margin
            if nwindow is not None:
                aligny += margin + (data.shape[0] - nwindow) // 2
                alignx += margin + (data.shape[1] - nwindow) // 2
            alignxlo = alignx + dx - nalignregion // 2
            alignxhi = alignx + dx + nalignregion // 2
            alignylo = aligny + dy - nalignregion // 2
            alignyhi = aligny + dy + nalignregion // 2

            aligndata = data[alignylo:alignyhi, alignxlo:alignxhi].copy()
            aligndata -= np.nanmedian(aligndata)
            aligndata = np.nan_to_num(aligndata, nan=0.0)
            max = np.unravel_index(np.argmax(aligndata, axis=None), aligndata.shape)
            ddy = max[0] - nalignregion // 2
            ddx = max[1] - nalignregion // 2
            print(
                "makeobject: maximum is offset by ddx = %+d px ddy = %+d px."
                % (ddx, ddy)
            )
            dx += ddx
            dy += ddy
            print("makeobject: refined offset is dx = %+d px dy = %+d px." % (dx, dy))
            if showalignment:
                tenue.image.show(aligndata, contrast=0.05)

        datashape = np.array(data.shape)
        newdata = np.full(datashape + 2 * margin, np.nan, dtype=float)
        xlo = margin - dx
        xhi = xlo + datashape[1]
        ylo = margin - dy
        yhi = ylo + datashape[0]
        newdata[ylo:yhi, xlo:xhi] = data
        data = newdata

        yc = int(data.shape[0] / 2)
        xc = int(data.shape[1] / 2)
        ys = int(yc - nwindow / 2)
        xs = int(xc - nwindow / 2)
        data = data[ys : ys + nwindow, xs : xs + nwindow]

        print("makeobject: subtracting sky.")
        data -= np.nanmedian(data, keepdims=True)

        return data

    print("makeobject: making %s object from %s." % (filter, directorypath))

    fitspathlist = tenue.path.getrawfitspaths(directorypath + "/object/", filter=filter)
    if len(fitspathlist) == 0:
        print("ERROR: no object files found.")
        return

    if refalpha == None or refdelta == None:
        print("makeobject: determining reference pointing.")
        pointinglist = list(readonepointing(fitspath) for fitspath in fitspathlist)
        refpointing = (np.max(pointinglist, axis=0) + np.min(pointinglist, axis=0)) / 2
        refalpha = refpointing[0]
        refdelta = refpointing[1]
    print(
        "makeobject: reference pointing is alpha = %.5f deg delta = %.5f deg."
        % (math.degrees(refalpha), math.degrees(refdelta))
    )

    if doskyimage:
        skystack = list(readonesky(fitspath) for fitspath in fitspathlist)
        print("makeobject: making sky image.")
        skymean, skysigma = tenue.image.clippedmeanandsigma(skystack, sigma=3, axis=0)
        sigma = tenue.image.clippedmean(skysigma, sigma=3) / math.sqrt(len(fitspathlist))
        print("makeobject: estimated noise in sky image is %.2f." % sigma)
        tenue.image.show(skymean, zscale=True)
        writesky(directorypath, skymean, filter, name="makeobject")
    else:
        skydata = 0

    objectstack = np.array(list(readoneobject(fitspath) for fitspath in fitspathlist))

    if sigma is None:
        print(
            "makeobject: averaging %d object files without rejection."
            % len(objectstack)
        )
        objectdata = np.average(objectstack, axis=0)
    else:
        print(
            "makeobject: averaging %d object files with rejection." % len(objectstack)
        )
        objectdata = tenue.image.clippedmean(objectstack, sigma=10, axis=0)

    tenue.image.show(objectdata, zscale=True, contrast=0.1)

    writeobject(directorypath, objectdata, filter, name="makeobject")

    print("makeobject: finished.")

    return
