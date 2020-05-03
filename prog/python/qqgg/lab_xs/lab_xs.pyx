from libc.math cimport tanh, atanh, sqrt

def diff_xs_eta(double e_proton, double charge, double eta, double x_1, double x_2):
    """Calculates the differential cross section as a function of the
    cosine of the pseudo rapidity eta of one photon in units of 1/GeV².
    Here dΩ=detadφ
    :param e_proton: proton energy per beam [GeV]
    :param charge: charge of the quark
    :param x_1: momentum fraction of the first quark
    :param x_2: momentum fraction of the second quark
    :param eta: pseudo rapidity
    :return: the differential cross section [GeV^{-2}]
    """
    cdef double rap = atanh((x_1 - x_2) / (x_1 + x_2))
    return (
        charge ** 4
        / (137.036 * e_proton) ** 2
        / (24 * x_1 * x_2)
        * ((tanh(eta - rap)) ** 2 + 1)
    )

def pT(double e_proton, double eta, double x_1, double x_2):
    cdef double tanh_eta = tanh(eta)
    return (
        2
        * e_proton
        * x_1
        * x_2
        / (x_1 + x_2 - (x_1 - x_2) * tanh_eta)
        * sqrt(1 - tanh_eta ** 2)
    )


def averaged_tchanel_q2(double e_proton, double x_1, double x_2):
    return 2 * x_1 * x_2 * e_proton ** 2

def second_eta(double eta, double x_1, double x_2):
    cdef double rap = atanh((x_1 - x_2) / (x_1 + x_2))
    return (-eta + 2 * rap)
