"""
Implementation of the analytical cross section for q q_bar ->
γγ in the lab frame.

Author: Valentin Boettcher <hiro@protagon.space>
"""

import numpy as np
import monte_carlo
import lhapdf
from numba import jit, vectorize, float64


@vectorize([float64(float64, float64, float64, float64)], nopython=True)
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
    x_1 = np.asarray(x_1)
    x_2 = np.asarray(x_2)
    cosθ = np.asarray(cosθ)
    assert (
        x_1.shape == x_2.shape == cosθ.shape
    ), "Invalid shapes for the event parameters."

    q_1 = (
        e_proton
        * x_1
        * np.array(
            [
                np.ones_like(cosθ),
                np.zeros_like(cosθ),
                np.zeros_like(cosθ),
                np.ones_like(cosθ),
            ]
        )
    )
    q_2 = (
        e_proton
        * x_2
        * np.array(
            [
                np.ones_like(cosθ),
                np.zeros_like(cosθ),
                np.zeros_like(cosθ),
                -np.ones_like(cosθ),
            ]
        )
    )
    g_3 = (
        2
        * e_proton
        * x_1
        * x_2
        / (2 * x_2 + (x_1 - x_2) * (1 - cosθ))
        * np.array(
            [1 * np.ones_like(cosθ), np.sqrt(1 - cosθ ** 2), np.zeros_like(cosθ), cosθ]
        )
    )
    g_4 = q_1 + q_2 - g_3

    q_1 = q_1.reshape(4, cosθ.size).T
    q_2 = q_2.reshape(4, cosθ.size).T
    g_3 = g_3.reshape(4, cosθ.size).T
    g_4 = g_4.reshape(4, cosθ.size).T

    return np.array([q_1, q_2, g_3, g_4])


@vectorize([float64(float64, float64, float64, float64, float64)], nopython=True)
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


@vectorize([float64(float64, float64, float64, float64, float64)], nopython=True)
def diff_xs_η(e_proton, charge, η, x_1, x_2):
    """Calculates the differential cross section as a function of the
    cosine of the pseudo rapidity η of one photon in units of 1/GeV².

    Here dΩ=dηdφ

    :param e_proton: proton energy per beam [GeV]
    :param charge: charge of the quark
    :param x_1: momentum fraction of the first quark
    :param x_2: momentum fraction of the second quark
    :param η: pseudo rapidity

    :return: the differential cross section [GeV^{-2}]
    """
    tanh_η = np.tanh(η)
    f = energy_factor(e_proton, charge, x_1, x_2)

    return (x_1 ** 2 * (1 - tanh_η) ** 2 + x_2 ** 2 * (1 + tanh_η) ** 2) / (
        x_1 * (1 - tanh_η) + x_2 * (1 + tanh_η)
    )


@vectorize([float64(float64, float64, float64)], nopython=True)
def averaged_tchanel_q2(e_proton, x_1, x_2):
    return 2 * x_1 * x_2 * e_proton ** 2

from numba.extending import get_cython_function_address

def get_xs_distribution_with_pdf(xs, q, e_hadron, quarks=None, pdf=None):
    """Creates a function that takes an event (type np.ndarray) of the
    form [cosθ, impulse fractions of quarks in hadron 1, impulse
    fractions of quarks in hadron 2] and returns the differential
    cross section for such an event. I would have used an object as
    argument, wasn't for the sampling function that needs a vector
    valued function. Cosθ can actually be any angular-like parameter
    as long as the xs has the corresponding parameter.

    :param xs: cross section function with signature (energy hadron, cosθ, x_1, x_2)
    :param q2: the momentum transfer Q^2 as a function with the signature
    (e_hadron, x_1, x_2)
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

    xfxQ2 = pdf.xfxQ2

    # @jit(float64(float64[4])) Unfortunately that does not work as yet!
    def distribution(event: np.ndarray) -> float:
        cosθ, x_1, x_2 = event

        q2_value = q(e_hadron, x_1, x_2)
        result = 0

        for quark, charge in quarks:
            xs_value = xs(e_hadron, charge, cosθ, x_1, x_2)
            result += (
                xfxQ2(quark, x_1, q2_value)
                / x_1
                * xfxQ2(quark, x_2, q2_value)
                / x_2
                * xs_value
            )

        return result

    return distribution, (pdf.xMin, pdf.xMax)

def sample_momenta(num_samples, dist, interval, e_hadron, upper_bound=None):
    res, eff = monte_carlo.sample_unweighted_array(
        num_samples, dist, interval, upper_bound=upper_bound, report_efficiency=True
    )
    cosθ, x_1, x_2 = res.T
    return momenta(e_hadron, x_1[None, :], x_2[None, :], cosθ[None, :]), eff
