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
    ax, points, dist, increment_borders, integral=None, color="orange", *args, **kwargs
):
    """Plot the distribution with VEGAS weights applied.

    :param ax: axis
    :param points: points
    :param dist: distribution
    :param increment_borders: increment borders
    """

    num_increments = increment_borders.size
    weighted_dist = dist.copy()

    var = 0
    total_weight = points.max() - points.min()
    for left_border, right_border in zip(increment_borders[:-1], increment_borders[1:]):
        length = right_border - left_border
        mask = (left_border <= points) & (points <= right_border)
        weighted_dist[mask] = dist[mask] * num_increments * length
        if integral:
            var += (
                np.sum((integral - weighted_dist[mask]) ** 2)
                / (weighted_dist[mask].size - 1)
                * length
                / total_weight
            )

    if integral:
        std = np.sqrt(var)
        ax.axhline(weighted_dist.mean(), color=color, linestyle="--")
        ax.axhspan(
            weighted_dist.mean() - std,
            weighted_dist.mean() + std,
            color=color,
            alpha=0.2,
            linestyle="--",
        )

    ax.plot(
        points, weighted_dist, *args, color=color, **kwargs,
    )


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

import matplotlib.gridspec as gridspec
from scipy.stats import chi2


def test_compatibility(hist_1, hist_2):
    hist_2.reshape(hist_1.shape)
    test = np.sum(
        np.divide(
            ((hist_1 - hist_2) ** 2).astype(float),
            (hist_1 + hist_2).astype(float),
            out=np.zeros_like(hist_1).astype(float),
            where=(hist_1 + hist_2) > 0,
        )
    )

    return test, 1 - chi2.cdf(test, len(hist_1))


def draw_ratio_plot(histograms, normalize_to=1, **kwargs):
    fig, (ax_hist, ax_ratio) = set_up_plot(
        2, 1, sharex=True, gridspec_kw=dict(height_ratios=[3, 1], hspace=0), **kwargs
    )

    reference, edges = histograms[0]["hist"]
    reference_error = np.sqrt(reference)

    ref_int = hist_integral(histograms[0]["hist"])
    reference = reference / ref_int
    reference_error = reference_error / ref_int

    for i, histogram in enumerate(histograms):
        heights, _ = (
            histogram["hist"]
            if "hist" in histogram
            else np.histogram(histogram["samples"], bins=edges)
        )

        if i > 0:
            test, ratio = test_compatibility(heights, histograms[0]["hist"][0])
            histogram["hist_kwargs"]["label"] = (
                histogram["hist_kwargs"].get("label", "")
                + "\n"
                + rf"$P = {ratio:.2f}\; T = {test:.2f}$"
            )

        integral = hist_integral([heights, edges])
        errors = np.sqrt(heights) / integral
        heights = heights / integral

        draw_histogram(
            ax_hist,
            [heights, edges],
            errorbars=errors,
            hist_kwargs=(
                histogram["hist_kwargs"] if "hist_kwargs" in histogram else dict()
            ),
            errorbar_kwargs=(
                histogram["errorbar_kwargs"]
                if "errorbar_kwargs" in histogram
                else dict()
            ),
            normalize_to=normalize_to,
        )

        set_up_axis(ax_ratio, pimp_top=False)
        draw_histogram(
            ax_ratio,
            [
                np.divide(
                    heights, reference, out=np.ones_like(heights), where=reference != 0
                ),
                edges,
            ],
            errorbars=np.divide(
                errors, reference, out=np.zeros_like(heights), where=reference != 0
            ),
            hist_kwargs=(
                histogram["hist_kwargs"] if "hist_kwargs" in histogram else dict()
            ),
            errorbar_kwargs=(
                histogram["errorbar_kwargs"]
                if "errorbar_kwargs" in histogram
                else dict()
            ),
            normalize_to=None,
        )

        ax_ratio.set_ylabel("ratio")

    return fig, (ax_hist, ax_ratio)


def hist_integral(hist):
    heights, edges = hist
    return heights @ (edges[1:] - edges[:-1])


def draw_histogram(
    ax,
    histogram,
    errorbars=True,
    hist_kwargs=dict(color="#1f77b4"),
    errorbar_kwargs=dict(),
    autoau=True,
    normalize_to=None,
):
    """Draws a histogram with optional errorbars using the step style.

    :param ax: axis to draw on
    :param histogram: an array of the form [heights, edges]
    :param hist_kwargs: keyword args to pass to `ax.step`
    :param errorbar_kwargs: keyword args to pass to `ax.errorbar`
    :param autoau: if set, the y axis will receive an a.u. label
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

    hist_plot = ax.step(edges, [heights[0], *heights], **hist_kwargs)

    if errorbars is not False:
        if "color" not in errorbar_kwargs:
            errorbar_kwargs["color"] = hist_plot[0].get_color()

        ax.errorbar(centers, heights, deviations, linestyle="none", **errorbar_kwargs)

    ax.set_xlim(*[edges[0], edges[-1]])
    ax.set_ylabel("a.u.")

    return ax


def draw_histo_auto(points, xlabel, bins=50, range=None, rethist=False, **kwargs):
    """Creates a histogram figure from sample points, normalized to unity.

    :param points: samples
    :param xlabel: label of the x axis
    :param bins: number of bins
    :param range: the range of the values
    :param rethist: whether to return the histogram as third argument
    :returns: figure, axis
    """

    hist = np.histogram(points, bins, range=range, **kwargs)
    fig, ax = set_up_plot()
    draw_histogram(ax, hist, normalize_to=1)

    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")

    return (fig, ax, hist) if rethist else (fig, ax)

def yoda_to_numpy(histo):
    edges = histo.xEdges()
    heights = np.array([bi.numEntries() for bi in histo])

    return heights, edges


def draw_yoda_histo_auto(h, xlabel, **kwargs):
    hist = yoda_to_numpy(h)
    fig, ax = set_up_plot()
    draw_histogram(ax, hist, errorbars=True, normalize_to=1, **kwargs)

    ax.set_xlabel(xlabel)
    return fig, ax
