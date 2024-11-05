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
    return mean.astype("float32")


def clippedsigma(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    return sigma.astype("float32")


def clippedmeanandsigma(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    return mean.astype("float32"), sigma.astype("float32")

def medianfilter(data, size):
     return scipy.ndimage.median_filter(data, size)

def uniformfilter(data, size):
    return scipy.ndimage.filters.uniform_filter(data, size=size, mode="nearest")


def show(data, zrange=False, zscale=False, contrast=0.25, zmin=None, zmax=None):

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

    tickinterval = int(math.pow(2, int(math.log2(np.max(data.shape) / 8))))
    ticks = list(range(0, np.max(data.shape), tickinterval))

    plt.figure(figsize=(10, 10))
    plt.imshow(data, origin="lower", norm=norm)
    plt.xticks(ticks)
    plt.yticks(ticks)
    plt.colorbar(fraction=0.046, pad=0.035)
    plt.show()
