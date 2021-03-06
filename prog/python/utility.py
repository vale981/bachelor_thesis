import matplotlib
import matplotlib.pyplot as plt
from scipy.constants import hbar, c, electron_volt
import matplotlib.ticker as ticker
import functools
import numpy as np
import os
import errno

###############################################################################
#                                   Utility                                   #
###############################################################################


def numpy_cache(cache_arg_name):
    """Unobstrusively creates a cache for a function returning a numpy
    array.  A keyword argument with the name `cache_arg_name` will be
    injected into the function. If that argument is not specified, the
    function will be evaluated without cache. Otherwise a wrapper will
    look for the cache file specified by that keyword argument
    (without extension) and loads the results from that file. If the
    file is not found, it will be generated by a call to the decorated
    function.

    Caution: No cache invalidation is performed!

    :param cache_arg_name: the argument name for the cache

    """

    def decorator(f):
        @functools.wraps(f)
        def caching_wrapper(*args, **kwargs):
            if cache_arg_name in kwargs:
                if not kwargs[cache_arg_name]:
                    del kwargs[cache_arg_name]
                    return f(*args, **kwargs)

                path = kwargs[cache_arg_name] + ".npy"
                if os.path.isfile(path):
                    name, result = np.load(path, allow_pickle=True)
                    print("Loading Cache: ", *name)
                    if f.__name__ == name[0]:
                        return result

                    raise RuntimeError(
                        "Trying to read to cache from another function:" + name
                    )

                del kwargs[cache_arg_name]

                result = f(*args, **kwargs)

                if not os.path.exists(os.path.dirname(path)):
                    try:
                        os.makedirs(os.path.dirname(path))
                    except OSError as exc:  # Guard against race condition
                        if exc.errno != errno.EEXIST:
                            raise

                np.save(path, ([f.__name__], result))

                return result

            return f(*args, **kwargs)

        return caching_wrapper

    return decorator


def minkowski_product(v_1, v_2):
    """Performs the standard 4-vector product between `v_1` and `v_2`"""
    assert v_1.shape == v_2.shape, "Invalid input shapes."

    if len(v_1.shape) == 1 or len(v_1.shape) == 1:
        v_1 = v_1.reshape(4)
        v_2 = v_2.reshape(4)

        return v_1[0] * v_2[0] - (v_1[1:] @ v_2[1:])

    prod = v_1[:, 0] * v_2[:, 0] - np.sum(v_1[:, 1:] * v_2[:, 1:], axis=1)

    if prod.size == 1:
        return prod.item(0)
    return prod


def gev_to_pb(xs):
    """Converts a cross section from 1/GeV^2 to pb."""
    return xs / (electron_volt ** 2) * (hbar * c) ** 2 * 1e22


def pb_to_gev(xs):
    """Converts a cross section from pb to 1/GeV^2."""
    return xs * (electron_volt ** 2) / ((hbar * c) ** 2 * 1e22)


def θ_to_η(θ):
    θ = np.asarray(θ)
    return -np.log(np.tan(θ / 2))


def η_to_θ(η):
    η = np.asarray(η)
    return 2 * np.arctan(np.exp(-η))


def η_to_pt(η, p):
    return p / np.cosh(η)


def tex_value(val, err=None, unit=None, prefix="", suffix="", prec=0, save=None):
    """Generates LaTeX output of a value with units and error."""

    if err:
        val, err, prec = scientific_round(val, err, retprec=True)
    else:
        val = np.round(val, prec)

    if prec == 0:
        val = int(val)
        if err:
            err = int(err)

    val_string = fr"{val:.{prec}f}" if prec > 0 else str(val)
    if err:
        val_string += fr"\pm {err:.{prec}f}" if prec > 0 else str(err)

    ret_string = r"\(" + prefix

    if unit is None:
        ret_string += val_string
    else:
        ret_string += fr"\SI{{{val_string}}}{{{unit}}}"

    ret_string += suffix + r"\)"

    if save is not None:
        os.makedirs(save[0], exist_ok=True)
        with open(f"{save[0]}/{save[1]}", "w") as f:
            f.write(ret_string)

    return ret_string


###############################################################################
#                                  Plot Porn                                  #
###############################################################################

matplotlib.rcParams.update(
    {
        "font.family": "serif",
        "text.usetex": False,
        "pgf.rcfonts": False,
        "lines.linewidth": 1,
    }
)


def pinmp_ticks(axis, ticks):
    axis.set_major_locator(ticker.MaxNLocator(ticks))
    axis.set_minor_locator(ticker.MaxNLocator(ticks * 10))
    return axis


def set_up_axis(ax, ticks=4, pimp_top=True):
    pinmp_ticks(ax.xaxis, ticks)
    pinmp_ticks(ax.yaxis, ticks)

    ax.grid(which="minor", alpha=0.1, linewidth=0.2)
    ax.grid(which="major", alpha=0.3, linewidth=0.2)

    ax.tick_params(right=True, top=pimp_top, which="both")

    return ax


@functools.wraps(plt.subplot)
def set_up_plot(*args, ticks=4, pimp_top=True, **kwargs):
    fig, axes = plt.subplots(*args, squeeze=True, **kwargs)

    if hasattr(axes, "__len__"):
        shape = axes.shape

        axes = np.array(
            [set_up_axis(ax, ticks, pimp_top) for ax in axes.flatten()]
        ).reshape(shape)

    else:
        axes = set_up_axis(axes)

    return fig, axes


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def save_fig(fig, title, folder="unsorted", size=(5, 4)):
    fig.set_size_inches(*size)
    fig.tight_layout()

    size = cm2inch(*size)
    os.makedirs(f"./figs/{folder}/", exist_ok=True)

    fig.savefig(f"./figs/{folder}/{title}.pdf")
    fig.savefig(f"./figs/{folder}/{title}.pgf")


def scientific_round(val, *err, retprec=False):
    """Scientifically rounds the values to the given errors."""
    val, err = np.asarray(val), np.asarray(err)
    if len(err.shape) == 1:
        err = np.array([err])
        err = err.T
    err = err.T

    if err.size == 1 and val.size > 1:
        err = np.ones_like(val) * err

    if len(err.shape) == 0:
        err = np.array([err])

    if val.size == 1 and err.shape[0] > 1:
        val = np.ones_like(err) * val

    i = np.floor(np.log10(err))
    first_digit = (err // 10 ** i).astype(int)
    prec = (-i + np.ones_like(err) * (first_digit <= 3)).astype(int)
    prec = np.max(prec, axis=1)

    def smart_round(value, precision):
        value = np.round(value, precision)
        if precision <= 0:
            value = value.astype(int)
        return value

    if val.size > 1:
        rounded = np.empty_like(val)
        rounded_err = np.empty_like(err)
        for n, (value, error, precision) in enumerate(zip(val, err, prec)):
            rounded[n] = smart_round(value, precision)
            rounded_err[n] = smart_round(error, precision)

        if retprec:
            return rounded, rounded_err, prec
        else:
            return rounded, rounded_err

    else:
        prec = prec[0]
        if retprec:
            return (smart_round(val, prec), *smart_round(err, prec)[0], prec)
        else:
            return (smart_round(val, prec), *smart_round(err, prec)[0])
