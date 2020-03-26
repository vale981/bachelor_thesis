"""
Implementation of the analytical cross section for q\bar{q} ->
\gamma\gamma

Author: Valentin Boettcher <hiro@protagon.space>
"""

import numpy as np
from scipy.constants import alpha

# NOTE: a more elegant solution would be a decorator
def energy_factor(charge, esp):
    """
    Calculates the factor common to all other values in this module

    Arguments:
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    return charge**4*(alpha/esp)**2/4


def diff_xs(theta, charge, esp):
    """
    Calculates the differential cross section as a function of the
    azimuth angle theta in units of 1/GeV^2.

    Arguments:
    theta -- azimuth angle
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f*(1 + 2/np.sin(theta)**2)

def diff_xs_eta(eta, charge, esp):
    """
    Calculates the differential cross section as a function of the
    pseudo rapidity of the photons in units of 1/GeV^2.

    Arguments:
    eta -- pseudo rapidity
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f*(1 + 2*np.cosh(eta)**2)

def total_xs_eta(eta, charge, esp):
    """
    Calculates the total cross section as a function of the pseudo
    rapidity of the photons in units of 1/GeV^2.  If the rapditiy is
    specified as a tuple, it is interpreted as an interval.  Otherwise
    the interval [-eta, eta] will be used.

    Arguments:
    eta -- pseudo rapidity (tuple or number)
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementar charge
    """

    f = energy_factor(charge, esp)
    if not isinstance(eta, tuple):
        eta = (-eta, eta)

    if len(eta) != 2:
        raise ValueError('Invalid eta cut.')

    def F(x):
        return -np.tanh(x) - np.log(1-np.tanh(x)) + np.log(1+np.tanh(x)) - x/4

    return -2*np.pi*f*(F(eta[0]) - F(eta[1]))
