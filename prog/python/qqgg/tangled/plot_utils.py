"""
Some shorthands for common plotting tasks related to the investigation
of monte-carlo methods in one rimension.

Author: Valentin Boettcher <hiro at protagon.space>
"""

import matplotlib.pyplot as plt
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
    deviations = np.sqrt(heights)

    if normalize_to is not None:
        integral = heights @ (edges[1:] - edges[:-1])
        heights = heights / integral * normalize_to
        deviations = deviations / integral * normalize_to

    ax.errorbar(centers, heights, deviations, linestyle="none", **errorbar_kwargs)
    ax.step(edges, [heights[0], *heights], **hist_kwargs)
    print([edges[0], edges[-1]])
    ax.set_xlim(*[edges[0], edges[-1]])

    return ax


def draw_histo_auto(points, xlabel, bins=50, range=None, **kwargs):
    """Creates a histogram figure from sample points, normalized to unity.

    :param points: samples
    :param xlabel: label of the x axis
    :param bins: number of bins
    :param range: the range of the values
    :returns: figure, axis
    """

    hist = np.histogram(points, bins, range=range, **kwargs)

    fig, ax = set_up_plot()
    draw_histogram(ax, hist, normalize_to=1)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")

    return fig, ax

def draw_yoda_histo(h, xlabel):
    edges = np.append(h.xMins(), h.xMax())
    heights = np.append(h.yVals(), h.yVals()[-1])
    centers = (edges[1:] + edges[:-1]) / 2

    fig, ax = set_up_plot()
    ax.errorbar(h.xVals(), h.yVals(), h.yErrs(), linestyle="none", color="orange")
    ax.step(edges, heights, color="#1f77b4", where="post")

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_xlim([h.xMin(), h.xMax()])
    return fig, ax
