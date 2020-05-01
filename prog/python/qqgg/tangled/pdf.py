"""
Implementation of the analytical cross section for q q_bar ->
γγ in the lab frame.

Author: Valentin Boettcher <hiro@protagon.space>
"""

import numpy as np
import monte_carlo
import lhapdf
from numba import jit, vectorize, float64, boolean
import lab_xs.lab_xs as c_xs


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


class Cut:
    def __init__(self):
        self._other = None
        self._current_comb = self._call

        self._greater_than = 0
        self._lower_than = np.inf

    def __gt__(self, greater_than):
        self._greater_than = greater_than

        return self

    def __lt__(self, lower_than):
        self._lower_than = lower_than

        return self

    def _or_comb(self, event):
        return self._call(event) or self._other(event)

    def _and_comb(self, event):
        return self._call(event) and self._other(event)

    def _call(self, event):
        return self._greater_than < self._calculate(event) < self._lower_than

    def _calculate(self, event):
        raise NotImplementedError('"_calulate" must be implemented.')

    def __call__(self, event):
        return self._current_comb(event)

    def __and__(self, other):
        self._other = other
        self._current_comb = self._and_comb

        return self

    def __or__(self, other):
        self._other = other
        self._current_comb = self._or_comb

        return self

    def apply(self, function):
        @wraps(function)
        def wrapper(event):
            if self(event):
                return function(event)

            return 0

        return wrapper


@vectorize([float64(float64, float64, float64)], nopython=True)
def averaged_tchanel_q2(e_proton, x_1, x_2):
    return 2 * x_1 * x_2 * e_proton ** 2


class CutpT(Cut):
    def __init__(self):
        super().__init__()

    def _calculate(self, event):
        e_hadron, eta, x_1, x_2 = event
        return c_xs.pT(e_hadron, eta, x_1, x_2)


class CutOtherEta(Cut):
    def __init__(self):
        super().__init__()

    def _calculate(self, event):
        _, η, x_1, x_2 = event
        rap = np.arctanh((x_1 - x_2) / (x_1 + x_2))
        return c_xs.second_eta(η, x_1, x_2)

def cached_pdf(pdf, q, points, e_hadron):
    x_min = pdf.xMin
    x_max = pdf.xMax
    Q2_max = 2 * e_hadron ** 2

    cache = np.array(
        [
            [
                pdf.xfxQ2(
                    q, xx := x_min + (x_max - x_min) * x / points, Q2_max / 100 * Q2
                )
                / xx
                for Q2 in range(100)
            ]
            for x in range(points)
        ]
    )

    def cached(x, q2):
        return cache[int((x - x_min) / (x_max - x_min) * points - 1)][
            int(q2 * 100 / Q2_max - 1)
        ]

    return cached


def get_xs_distribution_with_pdf(
    xs,
    q,
    e_hadron,
    quarks=None,
    pdf=None,
    cut=None,
    num_points_pdf=1000,
    vectorize=False,
):
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
    quarks = (
        quarks
        if quarks is not None
        else np.array(
            # [[5, -1 / 3], [4, 2 / 3], [3, -1 / 3], [2, 2 / 3], [1, -1 / 3]]
            [[1, -1 / 3]]
        )
    )  # all the light quarks

    supported_quarks = pdf.flavors()
    for flavor in quarks[:, 0]:
        assert flavor in supported_quarks, (
            "The PDF doesn't support the quark flavor " + flavor
        )

    xfxQ2 = pdf.xfxQ2

    def distribution(event: np.ndarray) -> float:
        if cut and not cut([e_hadron, *event]):
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

    def vectorized(events):
        result = np.empty(events.shape[0])
        for i in range(events.shape[0]):
            result[i] = distribution(events[i])
        return result

    return vectorized if vectorize else distribution, (pdf.xMin, pdf.xMax)

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
