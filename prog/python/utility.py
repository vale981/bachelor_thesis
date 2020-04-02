import matplotlib
import matplotlib.pyplot as plt
from SecondaryValue import SecondaryValue
from scipy.constants import hbar, c, electron_volt
import matplotlib.ticker as ticker
import numpy as np
import os
import tikzplotlib
###############################################################################
#                                   Utility                                   #
###############################################################################

def gev_to_pb(xs):
    """Converts a cross section from 1/GeV^2 to pb."""
    return xs/(electron_volt**2)*(hbar*c)**2*1e22

def θ_to_η(θ):
    θ = np.asarray(θ)
    return -np.log(np.tan(θ/2))

def η_to_θ(η):
    η = np.asarray(η)
    return 2*np.arctan(np.exp(-η))

def η_to_pt(η, p):
    return p/np.cosh(η)

def tex_value(val, err=None, unit='', prefix='', prec=10, save=None):
    """Generates LaTeX output of a value with units and error."""

    if err:
        val, err, prec = scientific_round(val, err, retprec=True)
    else:
        val = np.round(val, prec)

    val_string = fr'\({prefix}\SI{{{val:.{prec}f}'

    if err:
        val_string += fr'\pm {err:.{prec}f}'

    val_string += fr'}}{{{unit}}}\)'

    if save:
        os.makedirs(save[0], exist_ok=True)
        with open(f'{save[0]}/{save[1]}', 'w') as f:
            f.write(val_string)

    return val_string

###############################################################################
#                                  Plot Porn                                  #
###############################################################################

matplotlib.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,
    'pgf.rcfonts': False,
})

def pinmp_ticks(axis, ticks):
    axis.set_major_locator(ticker.MaxNLocator(ticks))
    axis.set_minor_locator(ticker.MaxNLocator(ticks*10))
    return axis

def set_up_plot(ticks=4, pimp_top=True, subplot=111, fig=None):
    if fig is None:
        fig = plt.figure()
    ax = fig.add_subplot(subplot)

    pinmp_ticks(ax.xaxis, ticks)
    pinmp_ticks(ax.yaxis, ticks)

    ax.grid(which='minor', alpha=.3)
    ax.grid(which='major', alpha=.5)


    if pimp_top:
        ax.tick_params(right=True, top=True, which='both')
    else:
        ax.tick_params(right=True, which='both')

    return fig, ax

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def save_fig(fig, title, folder='unsorted', size=(5, 4)):
    fig.set_size_inches(*size)
    fig.tight_layout()

    size = cm2inch(*size)
    os.makedirs(f'./figs/{folder}/', exist_ok=True)

    fig.savefig(f'./figs/{folder}/{title}.pdf')
    fig.savefig(f'./figs/{folder}/{title}.pgf')


def scientific_round(val, *err, retprec=False):
    """Scientifically rounds the values to the given errors."""
    val, err = np.asarray(val), np.asarray(err)
    if len(err.shape) == 1:
        err = np.array([err])
        err = err.T
    err = err.T

    if err.size == 1 and val.size > 1:
        err = np.ones_like(val)*err

    if len(err.shape) == 0:
        err = np.array([err])

    if val.size == 1 and err.shape[0] > 1:
        val = np.ones_like(err)*val

    i = np.floor(np.log10(err))
    first_digit = (err // 10**i).astype(int)
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
            return smart_round(val, prec), *smart_round(err, prec)[0], prec
        else:
            return smart_round(val, prec), *smart_round(err, prec)[0]
