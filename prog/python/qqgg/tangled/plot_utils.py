"""
Some shorthands for common plotting tasks related to the investigation
of monte-carlo methods in one rimension.

Author: Valentin Boettcher <hiro at protagon.space>
"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import yoda as yo
from collections.abc import Iterable
import yoda.plotting as yplt
import numpy as np
from utility import *


def plot_increments(ax, increment_borders, label=None, *args, **kwargs):
    """Plot the increment borders from a list.  The first and last one

    :param ax: the axis on which to draw
    :param list increment_borders: the borders of the increments
    :param str label: the label to apply to one of the vertical lines
    """

    ax.axvline(x=increment_borders[1], label=label, *args, **kwargs)

    for increment in increment_borders[1:-1]:
        ax.axvline(x=increment, *args, **kwargs)


def plot_vegas_weighted_distribution(
    ax, points, dist, increment_borders, *args, **kwargs
):
    """Plot the distribution with VEGAS weights applied.

    :param ax: axis
    :param points: points
    :param dist: distribution
    :param increment_borders: increment borders
    """

    num_increments = increment_borders.size
    weighted_dist = dist.copy()

    for left_border, right_border in zip(increment_borders[:-1], increment_borders[1:]):
        length = right_border - left_border
        mask = (left_border <= points) & (points <= right_border)
        weighted_dist[mask] = dist[mask] * num_increments * length

    ax.plot(points, weighted_dist, *args, **kwargs)


def plot_stratified_rho(ax, points, increment_borders, *args, **kwargs):
    """Plot the weighting distribution resulting from the increment
    borders.

    :param ax: axis
    :param points: points
    :param increment_borders: increment borders

    """

    num_increments = increment_borders.size
    ρ = np.empty_like(points)
    for left_border, right_border in zip(increment_borders[:-1], increment_borders[1:]):
        length = right_border - left_border
        mask = (left_border <= points) & (points <= right_border)
        ρ[mask] = 1 / (num_increments * length)

    ax.plot(points, ρ, *args, **kwargs)

def draw_histogram(
    ax,
    histogram,
    errorbars=True,
    hist_kwargs=dict(color="#1f77b4"),
    errorbar_kwargs=dict(color="orange"),
    normalize_to=None,
):
    """Draws a histogram with optional errorbars using the step style.

    :param ax: axis to draw on
    :param histogram: an array of the form [heights, edges]
    :param hist_kwargs: keyword args to pass to `ax.step`
    :param errorbar_kwargs: keyword args to pass to `ax.errorbar`
    :param normalize_to: if set, the histogram will be normalized to the value
    :returns: the given axis
    """

    heights, edges = histogram
    centers = (edges[1:] + edges[:-1]) / 2
    deviations = (
        (errorbars if isinstance(errorbars, (np.ndarray, list)) else np.sqrt(heights))
        if errorbars is not False
        else None
    )

    if normalize_to is not None:
        integral = hist_integral(histogram)
        heights = heights / integral * normalize_to
        if errorbars is not False:
            deviations = deviations / integral * normalize_to

    if errorbars is not False:
        ax.errorbar(centers, heights, deviations, linestyle="none", **errorbar_kwargs)

    ax.step(edges, [heights[0], *heights], **hist_kwargs)
    ax.set_xlim(*[edges[0], edges[-1]])

    return ax


@functools.wraps(yplt.plot)
def yoda_plot_wrapper(*args, **kwargs):
    fig, axs = yplt.plot(*args, **kwargs)
    if isinstance(axs, Iterable):
        axs = [set_up_axis(ax) for ax in axs if ax is not None]
    else:
        axs = set_up_axis(axs)

    return fig, axs

def samples_to_yoda(samples, bins, range=None, **kwargs):
    if range is None:
        range = [min(samples), max(samples)]

    hist = yo.Histo1D(bins, *range, **kwargs)
    for sample in samples:
        hist.fill(sample)

    return hist

def draw_histo_auto(
    points, xlabel, title="", bins=50, range=None, rethist=False, **kwargs
):
    """Creates a histogram figure from sample points, normalized to unity.

    :param points: samples
    :param xlabel: label of the x axis
    :param bins: number of bins
    :param range: the range of the values
    :param rethist: whether to return the histogram as third argument
    :returns: figure, axis
    """

    hist = yo.Histo1D(
        bins, *(range if range is not None else [min(points), max(points)]), title
    )

    for point in points:
        hist.fill(point)

    plot = yoda_plot_wrapper(hist, xlabel=xlabel, **kwargs)
    return (*plot, hist) if rethist else plot

def yoda_to_numpy(histo):
    histo.normalize(
        histo.numEntries() * ((histo.xMax() - histo.xMin()) / histo.numBins())
    )
    edges = np.append(histo.xMins(), histo.xMax())
    heights = histo.yVals().astype(int)

    return heights, edges


def draw_yoda_histo_auto(h, xlabel, **kwargs):
    hist = yoda_to_numpy(h)
    fig, ax = set_up_plot()
    draw_histogram(ax, hist, errorbars=True, normalize_to=1, **kwargs)

    ax.set_xlabel(xlabel)
    return fig, ax
