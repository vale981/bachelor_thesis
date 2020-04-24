import numpy as np
import matplotlib.pyplot as plt
import monte_carlo

"""
Implementation of the analytical cross section for q q_bar ->
gamma gamma

Author: Valentin Boettcher <hiro@protagon.space>
"""

import numpy as np

# NOTE: a more elegant solution would be a decorator
def energy_factor(charge, esp):
    """
    Calculates the factor common to all other values in this module

    Arguments:
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    return charge ** 4 / (137.036 * esp) ** 2 / 6


def diff_xs(θ, charge, esp):
    """
    Calculates the differential cross section as a function of the
    azimuth angle θ in units of 1/GeV².

    Here dΩ=sinθdθdφ

    Arguments:
    θ -- azimuth angle
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f * ((np.cos(θ) ** 2 + 1) / np.sin(θ) ** 2)


def diff_xs_cosθ(cosθ, charge, esp):
    """
    Calculates the differential cross section as a function of the
    cosine of the azimuth angle θ in units of 1/GeV².

    Here dΩ=d(cosθ)dφ

    Arguments:
    cosθ -- cosine of the azimuth angle
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f * ((cosθ ** 2 + 1) / (1 - cosθ ** 2))


def diff_xs_eta(η, charge, esp):
    """
    Calculates the differential cross section as a function of the
    pseudo rapidity of the photons in units of 1/GeV^2.

    This is actually the crossection dσ/(dφdη).

    Arguments:
    η -- pseudo rapidity
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    return f * (np.tanh(η) ** 2 + 1)


def diff_xs_p_t(p_t, charge, esp):
    """
    Calculates the differential cross section as a function of the
    transverse momentum (p_t) of the photons in units of 1/GeV^2.

    This is actually the crossection dσ/(dφdp_t).

    Arguments:
    p_t -- transverse momentum in GeV
    esp -- center of momentum energy in GeV
    charge -- charge of the particle in units of the elementary charge
    """

    f = energy_factor(charge, esp)
    sqrt_fact = np.sqrt(1 - (2 * p_t / esp) ** 2)
    return f / p_t * (1 / sqrt_fact + sqrt_fact)


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
        raise ValueError("Invalid η cut.")

    def F(x):
        return np.tanh(x) - 2 * x

    return 2 * np.pi * f * (F(η[0]) - F(η[1]))

@numpy_cache("momentum_cache")
def sample_momenta(sample_num, interval, charge, esp, seed=None, **kwargs):
    """Samples `sample_num` unweighted photon 4-momenta from the
    cross-section. Superflous kwargs are passed on to
    `sample_unweighted_array`.

    :param sample_num: number of samples to take
    :param interval: cosθ interval to sample from
    :param charge: the charge of the quark
    :param esp: center of mass energy
    :param seed: the seed for the rng, optional, default is system
        time

    :returns: an array of 4 photon momenta

    :rtype: np.ndarray

    """

    cosθ_sample = monte_carlo.sample_unweighted_array(
        sample_num, lambda x: diff_xs_cosθ(x, charge, esp), interval_cosθ, **kwargs
    )
    print(cosθ_sample)
    φ_sample = np.random.uniform(0, 1, sample_num)

    def make_momentum(esp, cosθ, φ):
        sinθ = np.sqrt(1 - cosθ ** 2)
        return np.array([1, sinθ * np.cos(φ), sinθ * np.sin(φ), cosθ],) * esp / 2

    momenta = np.array(
        [make_momentum(esp, cosθ, φ) for cosθ, φ in np.array([cosθ_sample, φ_sample]).T]
    )
    return momenta
