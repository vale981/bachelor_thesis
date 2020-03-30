"""
Implementation of the analytical cross section for q q_bar ->
gamma gamma

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

    return charge**4*(alpha/esp)**2/6


def diff_xs(θ, charge, esp):
    """
    Calculates the differential cross section as a function of the
    azimuth angle θ in units of 1/GeV².

    Arguments:
    θ -- azimuth angle
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f*((np.cos(θ)**2+1)/np.sin(θ)**2)

def diff_xs_cosθ(cosθ, charge, esp):
    """
    Calculates the differential cross section as a function of the
    cosine of the azimuth angle θ in units of 1/GeV².

    Arguments:
    θ -- azimuth angle
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f*((cosθ**2+1)/(1-cosθ**2))

def diff_xs_eta(η, charge, esp):
    """
    Calculates the differential cross section as a function of the
    pseudo rapidity of the photons in units of 1/GeV^2.

    Arguments:
    η -- pseudo rapidity
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f*(2*np.cosh(η)**2 - 1)

def total_xs_eta(η, charge, esp):
    """
    Calculates the total cross section as a function of the pseudo
    rapidity of the photons in units of 1/GeV^2.  If the rapditiy is
    specified as a tuple, it is interpreted as an interval.  Otherwise
    the interval [-η, η] will be used.

    Arguments:
    η -- pseudo rapidity (tuple or number)
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementar charge
    """

    f = energy_factor(charge, esp)
    if not isinstance(η, tuple):
        η = (-η, η)

    if len(η) != 2:
        raise ValueError('Invalid η cut.')

    def F(x):
        return np.tanh(x) - 2*x

    return 2*np.pi*f*(F(η[0]) - F(η[1]))
