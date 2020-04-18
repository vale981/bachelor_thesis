"""
Implementation of the analytical cross section for q q_bar ->
γγ in the lab frame.

Author: Valentin Boettcher <hiro@protagon.space>
"""

import numpy as np


def energy_factor(e_proton, charge, x_1, x_2):
    """Calculates the factor common to all other values in this module.

    :param e_proton: proton energy per beam
    :param charge: charge of the quark
    :param x_1: momentum fraction of the first quark
    :param x_2: momentum fraction of the second quark

    """
    return charge ** 4 / (137.036 * e_proton) ** 2 / (24 * x_1 * x_2)


def momenta(e_proton, x_1, x_2, cosθ):
    """Given the Energy of the incoming protons `e_proton` and the
    momentum fractions `x_1` and `x_2` as well as the cosine of the
    azimuth angle of the first photon the 4-momenta of all particles
    are calculated.
    """
    q_1 = e_proton * x_1 * np.array([1, 0, 0, 1])
    q_2 = e_proton * x_2 * np.array([1, 0, 0, -1])
    g_3 = (
        e_proton
        * x_1
        * x_2
        / (2 * x_2 + (x_1 - x_2) * (1 - cosθ))
        * np.array([1, cosθ, 0, sinθ])
    )
    g_4 = q_1 + q_2 - g_3

    return q_1, q_2, g_3, g_4


def diff_xs(e_proton, charge, cosθ, x_1, x_2):
    """Calculates the differential cross section as a function of the
    cosine of the azimuth angle θ of one photon in units of 1/GeV².

    Here dΩ=d(cosθ)dφ

    :param e_proton: proton energy per beam [GeV]
    :param charge: charge of the quark
    :param x_1: momentum fraction of the first quark
    :param x_2: momentum fraction of the second quark
    :param cosθ: the angle

    :return: the differential cross section [GeV^{-2}]
    """

    f = energy_factor(e_proton, charge, x_1, x_2)
    return (x_1 ** 2 * (1 - cosθ) ** 2 + x_2 ** 2 * (1 + cosθ) ** 2) / (
        (1 - cosθ ** 2) * (x_1 * (1 - cosθ) + x_2 * (1 + cosθ))
    )

def t_channel_q2(e_proton, cosθ, x_1, x_2):
    p, _, p_tag, _ = momenta(e_proton, x_1, x_2, cosθ)
    q = (p - p_tag)

    return -minkowski_product(q, q)
