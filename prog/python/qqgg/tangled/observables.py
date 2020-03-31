"""This module defines some observables on arrays of 4-pulses."""
import numpy as np

def p_t(p):
    """Transverse impulse

    :param p: array of 4-impulses
    """

    return np.linalg.norm(p[:,1:3], axis=1)

def η(p):
    """Pseudo rapidity.

    :param p: array of 4-impulses
    """

    return np.arccosh(np.linalg.norm(p[:,1:], axis=1)/p_t(p))*np.sign(p[:, 3])
