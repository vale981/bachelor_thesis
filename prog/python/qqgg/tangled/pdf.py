"""
Implementation of the analytical cross section for q q_bar ->
γγ in the lab frame.

Author: Valentin Boettcher <hiro@protagon.space>
"""

import numpy as np
import monte_carlo
import lhapdf

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
        * np.array([1, cosθ, 0, np.sqrt(1-cosθ**2)])
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

def get_xs_distribution_with_pdf(xs, q, e_hadron, quarks=None, pdf=None):
    """Creates a function that takes an event (type np.ndarray) of the
    form [cosθ, impulse fractions of quarks in hadron
    1, impulse fractions of quarks in hadron 2] and returns the
    differential cross section for such an event. I would have used an
    object as argument, wasn't for the sampling function that needs a
    vector valued function.

    :param xs: cross section function with signature (energy hadron, cosθ, x_1, x_2)
    :param q2: the momentum transfer Q^2 as a function with the signature
    (e_hadron, cosθ, x_1, x_2)
    :param quarks: the constituent quarks np.ndarray of the form [[id, charge], ...],
    the default is a proton
    :param pdf: the PDF to use, the default is "NNPDF31_lo_as_0118"
    :returns: differential cross section summed over flavors and weighted with the pdfs
    :rtype: function

    """

    pdf = pdf or lhapdf.mkPDF("NNPDF31_lo_as_0118", 0)
    quarks = quarks or np.array([[2, 2 / 3], [1, -1 / 3]])  # proton
    supported_quarks = pdf.flavors()
    for flavor in quarks[:, 0]:
        assert flavor in supported_quarks, (
            "The PDF doesn't support the quark flavor " + flavor
        )

    def distribution(event: np.ndarray) -> float:
        cosθ, x_1, x_2 = event

        q2_value = q(e_hadron, cosθ, x_1, x_2)
        result = 0

        for quark, charge in quarks:
            xs_value = xs(e_hadron, charge, cosθ, x_1, x_2)
            result += (
                pdf.xfxQ2(quark, x_1, q2_value)
                / x_1
                * pdf.xfxQ2(quark, x_2, q2_value)
                / x_2
                * xs_value
            )

        return result

    return distribution, (pdf.xMin, pdf.xMax)

def sample_momenta(num_samples, dist, interval, e_hadron):
    cosθ, x_1, x_2 = monte_carlo.sample_unweighted_array(
        num_samples, dist, interval
    ).T

    print(cosθ, x_1, x_2)

    return momenta(e_hadron, x_1, x_2, cosθ)[2:]  # only final state...
