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
    objectname,
    filter,
    path="{objectname}-{filter}.fits",
    starttimestamp=None,
    endtimestamp=None,
    exposuretime=None,
    gain=None,
    name="writeobject",
):
    path = path.format(objectname=objectname, filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(
        path,
        _objectdata,
        filter=filter,
        starttimestamp=starttimestamp,
        endtimestamp=endtimestamp,
        exposuretime=exposuretime,
        gain=gain,
    )
    return


def writesky(
    objectname,
    filter,
    path="{objectname}-{filter}-sky.fits",
    exposuretime=None,
    name="writesky",
):
    path = path.format(objectname=objectname, filter=filter)
    print("%s: writing %s." % (name, path))
    tenue.fits.writeproduct(path, _skydata, filter=filter, exposuretime=exposuretime)
    return


def _makesky(
    headerlist, datalist, filter=filter, skyclip=None, skystretch=None, objectname=None
):

    print("makeobject: making sky image.")

    newdatalist = []
    for data in datalist:
        newdata = np.copy(data)
        newdata -= np.nanmedian(newdata)
        if skyclip is not None:
            newdata[np.where(newdata >= +skyclip)] = np.nan
            newdata[np.where(newdata <= -skyclip)] = np.nan
        newdatalist.append(newdata)

    print("makeobject: averaging %d sky files with rejection." % len(datalist))
    global _skydata
    _skydata, skysigma = tenue.image.clippedmeanandsigma(newdatalist, sigma=3, axis=0)
    sigma = tenue.image.clippedmean(skysigma, sigma=3) / math.sqrt(len(datalist))
    totalexposuretime = np.sum(
        (tenue.instrument.exposuretime(header) for header in headerlist)
    )
    print(
        "makeobject: total exposure time in sky images is %.0f seconds."
        % totalexposuretime
    )
    print("makeobject: estimated noise in sky image is %.2f DN." % sigma)
    if skystretch is not None:
        tenue.image.show(
            _skydata / np.nanmedian(_skydata),
            zmin=1.0 - 0.5 * skystretch,
            zmax=1.0 + 0.5 * skystretch,
        )
    writesky(objectname, filter, exposuretime=totalexposuretime, name="makeobject")


def makeobject(
    fitspaths,
    objectname,
    filter,
    fitspathsslice=None,
    align=None,
    nalignregion=40,
    refalpha=None,
    refdelta=None,
    sigma=None,
    nwindow=None,
    showwindow=True,
    doskyimage=False,
    skyclip=None,
    skystretch=0.10,
    residualimageclip=None,
    triggertime=None,
    rejectfraction=0.25,
):

    ############################################################################

    print("makeobject: making %s object from %s." % (filter, fitspaths))

    fitspathlist = tenue.path.getrawfitspaths(
        fitspaths, filter=filter, fitspathsslice=fitspathsslice
    )

    if len(fitspathlist) == 0:
        print("ERROR: no object files found.")
        return

    print("makeobject: reading %d images." % len(fitspathlist))

    headerlist = []
    datalist = []
    for fitspath in fitspathlist:
        header, data = tenue.cook.cook(
            fitspath,
            name="makeobject",
            dooverscan=True,
            dotrim=True,
            dobias=True,
            dodark=True,
            doflat=True,
            dorotate=False,
        )
        headerlist.append(header)
        datalist.append(data)

    ############################################################################

    if refalpha == None or refdelta == None:
        print("makeobject: determining reference pointing.")
        alphalist = []
        deltalist = []
        for header, fitspath in zip(headerlist, fitspathlist):
            alpha = math.radians(tenue.instrument.alpha(header))
            delta = math.radians(tenue.instrument.delta(header))
            print(
                "makeobject: pointing for %s is alpha = %.5f deg delta = %.5f deg."
                % (os.path.basename(fitspath), math.degrees(alpha), math.degrees(delta))
            )
            alphalist.append(alpha)
            deltalist.append(delta)
            refalpha = (np.min(alphalist) + np.max(alphalist)) / 2
            refdelta = (np.min(deltalist) + np.max(deltalist)) / 2
    print(
        "makeobject: reference pointing is alpha = %.5f deg delta = %.5f deg."
        % (math.degrees(refalpha), math.degrees(refdelta))
    )

    ############################################################################

    if doskyimage:
        _makesky(
            headerlist,
            datalist,
            filter,
            skyclip=skyclip,
            skystretch=skystretch,
            objectname=objectname,
        )

    rotateddatalist = []
    for header, data in zip(headerlist, datalist):
        data -= np.nanmedian(data)
        if doskyimage:
            data -= _skydata
        rotateddata = tenue.instrument.dorotate(header, data)
        rotateddatalist.append(rotateddata)
    datalist = rotateddatalist

    ############################################################################

    dxlist = []
    dylist = []

    for header, fitspath in zip(headerlist, fitspathlist):

        alpha = math.radians(tenue.instrument.alpha(header))
        delta = math.radians(tenue.instrument.delta(header))
        pixelscale = math.radians(tenue.instrument.pixelscale(header))
        rotation = math.radians(tenue.instrument.rotation(header))

        dalpha = (alpha - refalpha) / pixelscale * math.cos(refdelta)
        ddelta = (delta - refdelta) / pixelscale
        dx = +int(np.round(dalpha * math.cos(rotation) - ddelta * math.sin(rotation)))
        dy = -int(np.round(dalpha * math.sin(rotation) + ddelta * math.cos(rotation)))
        print(
            "makeobject: %s: raw offset is dx = %+3d px dy = %+3d px."
            % (os.path.basename(fitspath), dx, dy)
        )

        dxlist.append(dx)
        dylist.append(dy)

    ############################################################################

    margin = 512

    # Extract the window data.

    windowdatalist = []

    if align is not None:

        for data, dx, dy in zip(datalist, dxlist, dylist):

            aligny = align[0] - margin
            alignx = align[1] - margin
            if nwindow is not None:
                aligny += margin + (data.shape[0] - nwindow) // 2
                alignx += margin + (data.shape[1] - nwindow) // 2
            alignxlo = alignx + dx - nalignregion // 2
            alignxhi = alignx + dx + nalignregion // 2
            alignylo = aligny + dy - nalignregion // 2
            alignyhi = aligny + dy + nalignregion // 2

            windowdata = data[alignylo:alignyhi, alignxlo:alignxhi].copy()
            windowdata -= np.nanmedian(windowdata)
            windowdata = np.nan_to_num(windowdata, nan=0.0)

            windowdatalist.append(windowdata)

    ############################################################################

    # Refine the offset

    if align is not None:

        newdxlist = []
        newdylist = []

        for windowdata, dx, dy in zip(windowdatalist, dxlist, dylist):

            filteredwindowdata = tenue.image.medianfilter(np.copy(windowdata), 3)

            imax = np.unravel_index(
                np.argmax(windowdata, axis=None), filteredwindowdata.shape
            )
            ddy = imax[0] - nalignregion // 2
            ddx = imax[1] - nalignregion // 2
            dx += ddx
            dy += ddy
            print("makeobject: refined offset is dx = %+d px dy = %+d px." % (dx, dy))

            newdxlist.append(dx)
            newdylist.append(dy)

        dxlist = newdxlist
        dylist = newdylist

    ############################################################################

    # Determine the merit.

    meritlist = []

    if align is not None:

        for windowdata in windowdatalist:

            merit = float(
                np.max(tenue.image.medianfilter(windowdata, 3))
                / tenue.image.clippedsigma(windowdata, sigma=3)
            )
            meritlist.append(merit)

    ############################################################################

    # Reject based on merit.

    if len(meritlist) > 0:

        noriginal = len(datalist)
        meritlimit = np.percentile(meritlist, 100 * rejectfraction)
        print(
            "makeobject: rejecting fraction %.2f of images with merit less than %.1f."
            % (rejectfraction, meritlimit)
        )

        x = [-np.inf] + sorted(meritlist)
        y = np.arange(0, noriginal + 1) / noriginal

        plt.figure()
        plt.step(x, y, where="post", color="C0")
        plt.ylim(0, 1)
        plt.xlim(left=0)
        plt.yticks(np.linspace(0, 1, 11))
        plt.axvline(meritlimit, color="C1")
        plt.axhline(rejectfraction, color="C1")
        plt.xlabel("Merit")
        plt.ylabel("Cumulative Fraction")
        plt.show()

        fitspathlist = list(
            fitspath
            for fitspath, merit in zip(fitspathlist, meritlist)
            if merit >= meritlimit
        )
        datalist = list(
            data for data, merit in zip(datalist, meritlist) if merit >= meritlimit
        )
        headerlist = list(
            header
            for header, merit in zip(headerlist, meritlist)
            if merit >= meritlimit
        )
        windowdatalist = list(
            windowdata
            for windowdata, merit in zip(windowdatalist, meritlist)
            if merit >= meritlimit
        )
        dxlist = list(dx for dx, merit in zip(dxlist, meritlist) if merit >= meritlimit)
        dylist = list(dy for dy, merit in zip(dylist, meritlist) if merit >= meritlimit)
        meritlist = list(merit for merit in meritlist if merit >= meritlimit)

        nfinal = len(datalist)
        print("makeobject: accepted %d images out of %d." % (nfinal, noriginal))
        print(
            "makeobject: rejected %d images out of %d."
            % ((noriginal - nfinal), noriginal)
        )

    ############################################################################

    # Show the data.

    if align is not None:

        for fitspath, merit, windowdata in zip(fitspathlist, meritlist, windowdatalist):

            print(
                "makeobject: %s: merit is %.1f." % (os.path.basename(fitspath), merit)
            )

            if showwindow:
                tenue.image.show(windowdata, contrast=0.05, small=True)

    ############################################################################

    # Mask residual image

    if residualimageclip is not None:

        print("makeobject: masking residual image.")

        lastdata = datalist[0] * 0
        for data in datalist:
            mask = np.where(lastdata >= residualimageclip)
            lastdata = np.copy(data)
            data[mask] = np.nan

        datalist = datalist[1:]
        headerlist = headerlist[1:]
        dxlist = dxlist[1:]
        dylist = dylist[1:]

    ############################################################################

    # Produce the aligned data.

    aligneddatalist = []

    for data, dx, dy in zip(datalist, dxlist, dylist):

        datashape = np.array(data.shape)

        aligneddata = np.full(datashape + 2 * margin, np.nan, dtype="float32")
        xlo = margin - dx
        xhi = xlo + datashape[1]
        ylo = margin - dy
        yhi = ylo + datashape[0]
        aligneddata[ylo:yhi, xlo:xhi] = data

        yc = int(aligneddata.shape[0] / 2)
        xc = int(aligneddata.shape[1] / 2)
        ys = int(yc - nwindow / 2)
        xs = int(xc - nwindow / 2)
        aligneddata = aligneddata[ys : ys + nwindow, xs : xs + nwindow]

        aligneddata -= np.nanmedian(aligneddata, keepdims=True)
        # tenue.image.show(aligneddata, zscale=True, contrast=0.5)

        aligneddatalist.append(aligneddata)

    ############################################################################

    global _objectdata
    if sigma is None:
        print(
            "makeobject: averaging %d object files without rejection."
            % len(aligneddatalist)
        )
        _objectdata = np.average(aligneddatalist, axis=0)
    else:
        print(
            "makeobject: averaging %d object files with rejection."
            % len(aligneddatalist)
        )
        _objectdata, objectsigma = tenue.image.clippedmeanandsigma(
            aligneddatalist, sigma=10, axis=0
        )
        sigma = tenue.image.clippedmean(objectsigma, sigma=3) / math.sqrt(
            len(aligneddatalist)
        )
        print("makeobject: estimated noise in object image is %.2f DN." % sigma)

    tenue.image.show(_objectdata, zscale=True, contrast=0.1)

    ############################################################################

    # Determine time properties of the stack.

    totalexposuretime = np.sum(
        (tenue.instrument.exposuretime(header) for header in headerlist)
    )
    print("makeobject: total exposure time is %.0f seconds." % totalexposuretime)

    gain = tenue.instrument.gain(headerlist[0]) * math.sqrt(len(aligneddatalist))
    print("makeobject: effective gain is %.2f e/DN." % gain)

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

    ############################################################################

    writeobject(
        objectname,
        filter,
        starttimestamp=starttimestamp,
        endtimestamp=endtimestamp,
        exposuretime=totalexposuretime,
        gain=gain,
        name="makeobject",
    )

    print("makeobject: finished.")

    return
