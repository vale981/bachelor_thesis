"""
Simple monte carlo integration implementation.
Author: Valentin Boettcher <hiro@protagon.space>
"""
import numpy as np
from scipy.optimize import minimize_scalar

def integrate(f, interval, point_density=1000, seed=None, **kwargs):
    """Monte-Carlo integrates the functin `f` in an interval.

    :param f: function of one variable, kwargs are passed to it
    :param tuple interval: a 2-tuple of numbers, specifiying the
        integration range

    :returns: the integration result

    :rtype: float
    """

    assert len(interval) == 2, 'An interval has two endpoints'

    a, b = interval
    if b < a:
        a, b = b, a

    if seed:
        np.random.seed(seed)

    interval_length = (b - a)
    num_points = int(interval_length * point_density)
    points = np.random.uniform(a, b, num_points)
    sample = f(points, **kwargs)
    integral = np.sum(sample)/num_points*interval_length
    deviation = np.std(sample)/np.sqrt(num_points - 1)*interval_length
    return integral, deviation


def sample(f, interval, point_density=1000,
         upper_bound=None, seed=None, **kwargs):
    assert len(interval) == 2, 'An interval has two endpoints'

    a, b = interval
    if b < a:
        a, b = b, a

    np.random.seed(seed)

    interval_length = (b - a)
    num_points_x = int(interval_length*point_density)

    if not upper_bound:
        upper_bound = minimize_scalar(lambda *args: -f(*args, **kwargs),
                                      bounds=interval, method='bounded')
        if upper_bound.success:
            upper_bound = -upper_bound.fun
        else:
            raise RuntimeError('Could not find an upper bound.')

    points = np.random.uniform([a, 0], [b, 1], [num_points_x, 2])
    sample_points = points[:, 0] \
        [np.where(f(points[:, 0]) > points[:, 1]*upper_bound)]
    return sample_points
