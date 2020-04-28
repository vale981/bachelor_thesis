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


def momenta(e_proton, x_1, x_2, cosθ, φ=None):
    """Given the Energy of the incoming protons `e_proton` and the
    momentum fractions `x_1` and `x_2` as well as the cosine of the
    azimuth angle of the first photon the 4-momenta of all particles
    are calculated.
    """
    x_1 = np.asarray(x_1)
    x_2 = np.asarray(x_2)
    cosθ = np.asarray(cosθ)

    if φ is None:
        φ = 0
        cosφ = np.ones_like(cosθ)
        sinφ = 0

    else:
        if φ == "rand":
            φ = np.random.uniform(0, 2 * np.pi, cosθ.shape)
        else:
            φ = np.asarray(φ)
        sinφ = np.sin(φ)
        cosφ = np.cos(φ)

    assert (
        x_1.shape == x_2.shape == cosθ.shape
    ), "Invalid shapes for the event parameters."

    sinθ = np.sqrt(1 - cosθ ** 2)

    ones = np.ones_like(cosθ)
    zeros = np.zeros_like(cosθ)

    q_1 = e_proton * x_1 * np.array([ones, zeros, zeros, ones,])
    q_2 = e_proton * x_2 * np.array([ones, zeros, zeros, -ones,])
    g_3 = (
        2
        * e_proton
        * x_1
        * x_2
        / (x_1 + x_2 - (x_1 - x_2) * cosθ)
        * np.array([1 * np.ones_like(cosθ), sinθ * sinφ, cosφ * sinθ, cosθ])
    )
    g_4 = q_1 + q_2 - g_3

    q_1 = q_1.reshape(4, cosθ.size).T
    q_2 = q_2.reshape(4, cosθ.size).T
    g_3 = g_3.reshape(4, cosθ.size).T
    g_4 = g_4.reshape(4, cosθ.size).T

    return np.array([q_1, q_2, g_3, g_4])


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

    rap = np.arctanh((x_1 - x_2) / (x_1 + x_2))
    f = energy_factor(e_proton, charge, x_1, x_2)

    return f * ((np.tanh(η - rap)) ** 2 + 1)


@vectorize([float64(float64, float64, float64)], nopython=True)
def averaged_tchanel_q2(e_proton, x_1, x_2):
    return 2 * x_1 * x_2 * e_proton ** 2


def cut_pT_from_eta(greater_than=0):
    def cut(e_proton, η, x1, x2):
        cosθ = np.cos(η_to_θ(η))
        _, _, p1, p2 = momenta(e_proton, x1, x2, cosθ)
        return (
            np.sqrt((p1[0][1:3] ** 2).sum()) > greater_than
            and np.sqrt((p2[0][1:3] ** 2).sum()) > greater_than
        )

    return cut

from numba.extending import get_cython_function_address


def get_xs_distribution_with_pdf(xs, q, e_hadron, quarks=None, pdf=None, cut=None):
    """Creates a function that takes an event (type np.ndarray) of the
    form [angle_arg, impulse fractions of quarks in hadron 1, impulse
    fractions of quarks in hadron 2] and returns the differential
    cross section for such an event. I would have used an object as
    argument, wasn't for the sampling function that needs a vector
    valued function. Angle_Arg can actually be any angular-like parameter
    as long as the xs has the corresponding parameter.

    :param xs: cross section function with signature (energy hadron, angle_arg, x_1, x_2)
    :param q2: the momentum transfer Q^2 as a function with the signature
    (e_hadron, x_1, x_2)
    :param quarks: the constituent quarks np.ndarray of the form [[id, charge], ...],
    the default is a proton
    :param pdf: the PDF to use, the default is "NNPDF31_lo_as_0118"
    :param cut: cut function with signature (energy hadron, angle_arg, x_1,
    x_2) to return 0, when the event does not fit the cut

    :returns: differential cross section summed over flavors and weighted with the pdfs
    :rtype: function
    """

    pdf = pdf or lhapdf.mkPDF("NNPDF31_lo_as_0118", 0)
    quarks = quarks or np.array(
        [[5, -1 / 3], [4, 2 / 3], [3, -1 / 3], [2, 2 / 3], [1, -1 / 3]]
    )  # proton
    supported_quarks = pdf.flavors()
    for flavor in quarks[:, 0]:
        assert flavor in supported_quarks, (
            "The PDF doesn't support the quark flavor " + flavor
        )

    xfxQ2 = pdf.xfxQ2

    # @jit(float64(float64[4])) Unfortunately that does not work as yet!
    def distribution(event: np.ndarray) -> float:
        if cut and not cut(e_hadron, *event):
            return 0

        angle_arg, x_1, x_2 = event

        q2_value = q(e_hadron, x_1, x_2)
        result = 0

        for quark, charge in quarks:
            xs_value = xs(e_hadron, charge, angle_arg, x_1, x_2)
            result += (
                (xfxQ2(quark, x_1, q2_value) + xfxQ2(-quark, x_1, q2_value))
                / x_1
                * (xfxQ2(-quark, x_2, q2_value) + xfxQ2(quark, x_2, q2_value))
                / x_2
                * xs_value
            )

        return result

    return distribution, (pdf.xMin, pdf.xMax)

def sample_momenta(num_samples, dist, interval, e_hadron, upper_bound=None, **kwargs):
    res, eff = monte_carlo.sample_unweighted_array(
        num_samples,
        dist,
        interval,
        upper_bound=upper_bound,
        report_efficiency=True,
        **kwargs
    )
    cosθ, x_1, x_2 = res.T
    return momenta(e_hadron, x_1[None, :], x_2[None, :], cosθ[None, :]), eff
