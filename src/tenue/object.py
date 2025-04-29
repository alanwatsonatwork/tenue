import math
from datetime import datetime
import os.path
import numpy as np
import matplotlib.pyplot as plt
import warnings

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

    nz = len(datalist)
    ny = datalist[0].shape[0]
    nx = datalist[0].shape[1]
    newdatastack = np.full([nz, ny, nx], np.nan, dtype="float32")
    for iz in range(nz - 1):
        newdata = datalist[iz].copy()
        newdata -= np.nanmedian(newdata)
        if skyclip is not None:
            newdata[np.where(newdata >= +skyclip)] = np.nan
            newdata[np.where(newdata <= -skyclip)] = np.nan
        newdatastack[iz, :, :] = newdata

    print("makeobject: averaging %d sky files with rejection." % len(datalist))
    global _skydata
    _skydata, skysigma = tenue.image.clippedmeanandsigma(newdatastack, sigma=3, axis=0)
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
    objectname,
    fitspaths,
    filter,
    fitspathsslice=None,
    align=None,
    nalignmentwindow=40,
    refalpha=None,
    refdelta=None,
    sigma=None,
    ninputwindow=None,
    noutputwindow=None,
    showalignment=True,
    doskyimage=False,
    skyclip=None,
    skystretch=0.10,
    residualimageclip=None,
    triggertime=None,
    rejectfraction=0.25,
):
    """
    Make an object image by cooking, sky-subtracting, aligning, and combining a
    stack of raw images.

    :param objectname:
        The ``objectname`` argument must be a string. The output FITS files will
        be named ``objectname`` followed by ``"-"`` followed by ``filter``
        followed by one of ``".fits"`` or ``"-sky.fits"``.

    :param fitspaths:
        The ``fitspath``, ``filter``, and ``fitspathsslice`` arguments are
        passed to :func:`~tenue.path.getrawfitspaths` to generate the list of
        raw FITS files.

    :param filter:
        The ``filter`` argument is passed to :func:`~tenue.path.getrawfitspaths`
        to allow selecting raw FITS files in a given filter.

    :param fitspathsslice:
        The ``fitspathsslive``argument can be either a slice, either of the
        string ``"firsthalf"`` or ``"secondhalf"``, or ``None``. It is passed to
        :func:`~tenue.path.getrawfitspaths` to allow selecting a subset of raw
        FITS files.

    :param align:
        The ``align`` argument may be a list of two numbers or ``None``. If it
        is a list of two numbers, they define the position in the output image
        of the star that will be used for alignment refinement. If it is
        ``None``, the alignment will not be refined.

    :param nalignmentwindow:
        The ``nalignmentwindow`` argument may be a positive number. It defines
        the size of the window used to refine the alignment.

    :param refalpha:
    :param refdelta:
        The ``refalpha`` and ``refdelta`` may be the reference position for
        alignment in degrees or may be both ``None``. If they are not, the
        reference position is chosen as the average of the extremes of the
        unrefined pointings.

    :param sigma:
        The ``sigma`` argument may be a positive number or ``None``. If it is a
        positive number, the raw files will be combined with sigma-clipping at
        this level. If it is ``None``, they will be combined without
        sigma-clipping.

    :param ninputwindow:
        The ``ninputwindow`` argument may be a positive integer or ``None``. If
        it is a positive integer, the raw files are clipped to a window of this
        size centered on the data region. If it is ``None``, they are not.

    :param noutputwindow:
        The ``noutputwindow`` argument may be a positive integer or ``None``. If
        it is a positive integer, the output object file is clipped or expanded
        to a window of this size centered on the reference position. If it is
        ``None``, they are not.

    :param showalignment:
        The ``showalignment`` argument may be ``True`` or ``False``. If
        ``True``, the alignment window is shown for each non-rejected image.

    :param doskyimage:
        The ``doskyimage`` argument may be ``True`` or ``False``. If ``True``, a
        sky image is formed and subtracted.

    :param skyclip:
        The ``skyclip`` argument may be a positive number or ``None``. If it is
        a positive number, the individual sky files are clipped at this level
        before combining. If it is ``None``, they will be combined without
        clipping.

    :param skystretch:
        The stretch used to show the sky image.

    :param residualimageclip:
        The ``doskyimage`` argument may be a positive number or ``None``. If it
        is a positive number, then pixels at least this bright in the previous
        image (after cooking) are masked in the current. If it is ``None``, they
        are not.

    :param triggertime:
        The ``triggertime`` argument may be the trigger time as an ISO 8901 time
        represented as a string or ``None``. If it is not ``None``, it is used
        to calculate and print the start and end times of the observation
        relative to the trigger time.

    :param rejectfraction:
        The ``rejectfraction`` argument may be a number. When refining the
        alignment, each image is assigned a merit value equal to the maximum
        value in the alignment window divided by an estimator of the background
        noise in the alignment window. The images with the lowest merit, up to a
        fraction given by ``rejectfraction``, are rejected and not combined.

    :return: None
    """

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
            nwindow=ninputwindow,
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
    del datalist
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
        dx = +np.round(dalpha * math.cos(rotation) - ddelta * math.sin(rotation))
        dy = -np.round(dalpha * math.sin(rotation) + ddelta * math.cos(rotation))

        dx += tenue.instrument.boresightdx(header)
        dy += tenue.instrument.boresightdy(header)

        dx = int(dx + 0.5)
        dy = int(dy + 0.5)

        print(
            "makeobject: %s: raw offset is dx = %+3d px dy = %+3d px."
            % (os.path.basename(fitspath), dx, dy)
        )

        dxlist.append(dx)
        dylist.append(dy)

    ############################################################################

    # Extract the alignment window data.

    alignmentwindowdatalist = []

    if align is not None:

        for data, dx, dy in zip(datalist, dxlist, dylist):

            nydata = data.shape[0]
            nxdata = data.shape[1]

            aligny = align[0]
            alignx = align[1]
            alignylo = aligny + dy - nalignmentwindow // 2 + nydata // 2
            alignyhi = aligny + dy + nalignmentwindow // 2 + nydata // 2
            alignxlo = alignx + dx - nalignmentwindow // 2 + nxdata // 2
            alignxhi = alignx + dx + nalignmentwindow // 2 + nxdata // 2

            alignmentwindowdata = data[alignylo:alignyhi, alignxlo:alignxhi].copy()
            alignmentwindowdata -= np.nanmedian(windowdata)
            alignmentwindowdata = np.nan_to_num(windowdata, nan=0.0)

            alignmentwindowdatalist.append(alignmentwindowdata)

    ############################################################################

    # Refine the offset

    if align is not None:

        newdxlist = []
        newdylist = []

        for windowdata, dx, dy in zip(alignmentwindowdatalist, dxlist, dylist):

            filteredwindowdata = tenue.image.medianfilter(np.copy(windowdata), 3)

            imax = np.unravel_index(
                np.argmax(windowdata, axis=None), filteredwindowdata.shape
            )
            ddy = imax[0] - nalignmentwindow // 2
            ddx = imax[1] - nalignmentwindow // 2
            dx += ddx
            dy += ddy
            print("makeobject: refined offset is dx = %+d px dy = %+d px." % (dx, dy))

            newdxlist.append(dx)
            newdylist.append(dy)

        dxlist = newdxlist
        dylist = newdylist

    ############################################################################

    nmargin = max(np.max(np.abs(np.array(dxlist))), np.max(np.abs(np.array(dylist))))
    if noutputwindow is not None:
        nmargin = max(nmargin, int(noutputwindow - min(data.shape[0], data.shape[1])))
    print("makeobject: nmargin is %d." % nmargin)

    ############################################################################

    # Determine the merit.

    meritlist = []

    if align is not None:

        for windowdata in alignmentwindowdatalist:

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
            for windowdata, merit in zip(alignmentwindowdatalist, meritlist)
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

            if showalignment:
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

    # Produce the aligned and windowed data.

    aligneddatalist = []

    for data, dx, dy in zip(datalist, dxlist, dylist):

        nydata = data.shape[0]
        nxdata = data.shape[1]

        aligneddata = np.empty((nydata + 2 * nmargin, nxdata + 2 * nmargin))
        aligneddata.fill(np.nan)

        ylo = nmargin - dy
        yhi = ylo + nydata
        xlo = nmargin - dx
        xhi = xlo + nxdata

        aligneddata[ylo:yhi, xlo:xhi] = data

        if noutputwindow is not None:
            ylo = (aligneddata.shape[1] - noutputwindow) // 2
            yhi = ylo + noutputwindow
            xlo = (aligneddata.shape[0] - noutputwindow) // 2
            xhi = xlo + noutputwindow
            aligneddata = np.copy(aligneddata[ylo:yhi, xlo:xhi])

        aligneddata -= np.nanmedian(aligneddata, keepdims=True)

        aligneddatalist.append(aligneddata)

    ############################################################################

    global _objectdata
    if sigma is None:
        print(
            "makeobject: averaging %d object files without rejection."
            % len(aligneddatalist)
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            _objectdata = np.nanmean(aligneddatalist, axis=0)
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
        triggertimestamp = datetime.fromisoformat(triggertime + "+00:00").timestamp()
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
