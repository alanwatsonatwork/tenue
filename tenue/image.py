import math
import warnings

import numpy as np

import astropy.stats
import astropy.visualization

import scipy.ndimage

import matplotlib.pyplot as plt


def clippedmean(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    if isinstance(mean, np.ndarray):
        mean = mean.astype("float32")
    return mean


def clippedsigma(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    if isinstance(sigma, np.ndarray):
        sigma = sigma.astype("float32")
    return sigma


def clippedmeanandsigma(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    if isinstance(mean, np.ndarray):
        mean = mean.astype("float32")
    if isinstance(sigma, np.ndarray):
        sigma = sigma.astype("float32")
    return mean, sigma


def medianfilter(data, size):
    return scipy.ndimage.median_filter(data, size)


def uniformfilter(data, size):
    return scipy.ndimage.filters.uniform_filter(data, size=size, mode="nearest")


def show(
    data, zrange=False, zscale=False, contrast=0.25, zmin=None, zmax=None, small=False
):

    if zmin is not None and zmax is not None:
        interval = astropy.visualization.ManualInterval(zmin, zmax)
    elif zrange:
        interval = astropy.visualization.MinMaxInterval()
    else:
        interval = astropy.visualization.ZScaleInterval(contrast=contrast)
    stretch = astropy.visualization.LinearStretch()
    norm = astropy.visualization.ImageNormalize(
        data, interval=interval, stretch=stretch
    )

    if np.max(data.shape) > 1000:
        tickinterval = 100
    else:
        tickinterval = int(math.pow(2, int(math.log2(np.max(data.shape) / 16))))
    ticks = list(range(0, np.max(data.shape), tickinterval))

    if small:
        plt.figure(figsize=(5, 5))
    else:
        plt.figure(figsize=(10, 10))
    plt.imshow(data, origin="lower", norm=norm)
    plt.xticks(ticks, rotation=90)
    plt.yticks(ticks)
    plt.colorbar(fraction=0.046, pad=0.035)
    plt.show()
