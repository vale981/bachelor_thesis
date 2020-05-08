"""This module defines some observables on arrays of 4-pulses."""
import numpy as np
from utility import minkowski_product

def p_t(p):
    """Transverse momentum

    :param p: array of 4-momenta
    """

    return np.linalg.norm(p[:,1:3], axis=1)

def η(p):
    """Pseudo rapidity.

    :param p: array of 4-momenta
    """

    return np.arccosh(np.linalg.norm(p[:,1:], axis=1)/p_t(p))*np.sign(p[:, 3])

def inv_m(p_1, p_2):
    """Invariant mass off the final state system.

    :param p_1: array of 4-momenta, first fs particle
    :param p_2: array of 4-momenta, second fs particle
    """

    total_p = p_1 + p_2
    return np.sqrt(minkowski_product(total_p, total_p))

def cosθ(p):
    return p[:, 3] / p[:, 0]
