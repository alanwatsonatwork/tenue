import numpy as np

import astropy.stats

import warnings


def clippedmean(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    return mean


def clippedsigma(stack, sigma=3.0, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", Warning)
        mean, median, sigma = astropy.stats.sigma_clipped_stats(
            np.array(stack), sigma=sigma, axis=axis, cenfunc="median", stdfunc="mad_std"
        )
    return sigma
