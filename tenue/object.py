import math
from datetime import datetime
import os.path
import numpy as np
import matplotlib.pyplot as plt
import statistics

import tenue.cook
import tenue.fits
import tenue.image
import tenue.instrument
import tenue.path

_skydata = None
_objectdata = None


def writeobject(
    objectname, filter, path="{objectname}-{filter}.fits", name="writeobject"
):
    path = path.format(objectname=objectname, filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _objectdata, filter=filter)
    return


def writesky(
    objectname, filter, path="{objectname}-{filter}-sky.fits", name="writesky"
):
    path = path.format(objectname=objectname, filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _skydata, filter=filter)
    return


def makeobject(
    fitspaths,
    objectname,
    filter,
    align=None,
    nalignregion=40,
    refalpha=None,
    refdelta=None,
    sigma=None,
    nwindow=None,
    showalignment=True,
    doskyimage=False,
    skyclip=None,
    triggertime=None,
    rejectfraction=0.25,
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
        if skyclip is not None:
            data[np.where(data >= +skyclip)] = np.nan
            data[np.where(data <= -skyclip)] = np.nan
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
        if align is not None:

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

            merit = float(
                np.max(tenue.image.medianfilter(aligndata, 3))
                / tenue.image.clippedsigma(aligndata, sigma=3)
            )
            print("makeobject: merit is %.1f." % merit)
            meritlist.append(merit)

            if showalignment:
                tenue.image.show(aligndata, contrast=0.05)

        datashape = np.array(data.shape)
        newdata = np.full(datashape + 2 * margin, np.nan, dtype="float32")
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

    print("makeobject: making %s object from %s." % (filter, fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(fitspaths, filter=filter)
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

    global _skydata
    if doskyimage:
        skystack = list(readonesky(fitspath) for fitspath in fitspathlist)
        print("makeobject: making sky image.")
        _skydata, skysigma = tenue.image.clippedmeanandsigma(skystack, sigma=3, axis=0)
        sigma = tenue.image.clippedmean(skysigma, sigma=3) / math.sqrt(
            len(fitspathlist)
        )
        print("makeobject: estimated noise in sky image is %.2f DN." % sigma)
        _skydata[np.where(np.isnan(_skydata))] = 0
        tenue.image.show(_skydata, zmin=-20, zmax=50)
        writesky(objectname, filter, name="makeobject")
    else:
        _skydata = 0

    meritlist = []
    objectlist = list(readoneobject(fitspath) for fitspath in fitspathlist)
    headerlist = list(tenue.fits.readrawheader(fitspath) for fitspath in fitspathlist)
    
    if len(meritlist) > 0:

        plt.figure()
        n, bins = np.histogram(meritlist, bins=50)
        f = np.cumsum(n) / np.sum(n)
        plt.plot(bins[:-1], f)
        plt.ylim(0, 1)
        plt.xlim(left=0)
        plt.axhline(rejectfraction)
        plt.xlabel("Merit")
        plt.ylabel("Cumulative Fraction")
        plt.show()    

        noriginal = len(objectlist)
        meritlimit = np.percentile(meritlist, 100 * rejectfraction)
        print(
            "makeobject: rejecting fraction %.2f of images with merit less than %.1f."
            % (rejectfraction, meritlimit)
        )
        objectlist = list(
            object
            for object, merit in zip(objectlist, meritlist)
            if merit >= meritlimit
        )
        headerlist = list(
            header
            for header, merit in zip(headerlist, meritlist)
            if merit >= meritlimit
        )
        nfinal = len(objectlist)
        print("makeobject: accepted %d images out of %d." % (nfinal, noriginal))
        print("makeobject: rejected %d images out of %d." % ((noriginal - nfinal), noriginal))

    global _objectdata
    if sigma is None:
        print(
            "makeobject: averaging %d object files without rejection." % len(objectlist)
        )
        _objectdata = np.average(objectlist, axis=0)
    else:
        print("makeobject: averaging %d object files with rejection." % len(objectlist))
        _objectdata, objectsigma = tenue.image.clippedmeanandsigma(
            objectlist, sigma=10, axis=0
        )
        sigma = tenue.image.clippedmean(objectsigma, sigma=3) / math.sqrt(
            len(fitspathlist)
        )
        print("makeobject: estimated noise in object image is %.2f DN." % sigma)

    tenue.image.show(_objectdata, zscale=True, contrast=0.1)

    writeobject(objectname, filter, name="makeobject")

    # Determine time properties of the stack.

    starttimestamp = min(
        (tenue.instrument.starttimestamp(header) for header in headerlist)
    )
    print(
        "makeobject: observation start time is %s UTC."
        % datetime.utcfromtimestamp(starttimestamp).isoformat(" ", "seconds")
    )

    endtimestamp = max((tenue.instrument.endtimestamp(header) for header in headerlist))
    print(
        "makeobject: observation end   time is %s UTC."
        % datetime.utcfromtimestamp(endtimestamp).isoformat(" ", "seconds")
    )

    if triggertime is not None:
        triggertimestamp = datetime.fromisoformat(triggertime + "Z").timestamp()
        print(
            "makeobject: observations are from %.0f to %.0f seconds after the trigger"
            % ((starttimestamp - triggertimestamp), (endtimestamp - triggertimestamp))
        )
        print(
            "makeobject: observations are from %.2f to %.2f minutes after the trigger"
            % (
                (starttimestamp - triggertimestamp) / 60,
                (endtimestamp - triggertimestamp) / 60,
            )
        )
        print(
            "makeobject: observations are from %.2f to %.2f hours after the trigger"
            % (
                (starttimestamp - triggertimestamp) / 3600,
                (endtimestamp - triggertimestamp) / 3600,
            )
        )
        print(
            "makeobject: observations are from %.2f to %.2f days after the trigger"
            % (
                (starttimestamp - triggertimestamp) / 86400,
                (endtimestamp - triggertimestamp) / 86400,
            )
        )

    print("makeobject: finished.")

    return
