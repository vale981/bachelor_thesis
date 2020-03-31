"""
Simple monte carlo integration implementation.
Author: Valentin Boettcher <hiro@protagon.space>
"""
import numpy as np
from scipy.optimize import minimize_scalar
from functools import wraps


def _process_interval(interval):
    assert len(interval) == 2, 'An interval has two endpoints'

    a, b = interval
    if b < a:
        a, b = b, a

    return interval


@process_arg(1, _process_interval)
def integrate(f, interval, point_density=1000, seed=None, **kwargs):
    """Monte-Carlo integrates the functin `f` in an interval.

    :param f: function of one variable, kwargs are passed to it
    :param tuple interval: a 2-tuple of numbers, specifiying the
        integration range

    :returns: the integration result

    :rtype: float
    """

    interval = _process_interval(interval)
    if seed:
        np.random.seed(seed)

    interval_length = (interval[1] - interval[0])
    num_points = int(interval_length * point_density)

    points = np.random.uniform(interval[0], interval[1], num_points)

    sample = f(points, **kwargs)
    integral = np.sum(sample)/num_points*interval_length
    deviation = np.std(sample)/np.sqrt(num_points - 1)*interval_length

    return integral, deviation

def find_upper_bound(f, interval, **kwargs):
    """Find the upper bound of a function.

    :param f: function of one scalar and some kwargs that are passed
              on to it
    :param interval: interval to look in

    :returns: the upper bound of the function
    :rtype: float
    """

    upper_bound = minimize_scalar(lambda *args: -f(*args, **kwargs),
                                  bounds=interval, method='bounded')
    if upper_bound.success:
        return -upper_bound.fun
    else:
        raise RuntimeError('Could not find an upper bound.')

def sample_unweighted(f, interval, upper_bound=None, seed=None,
                    chunk_size=100, **kwargs):
    """Samples a distribution proportional to f by hit and miss.
    Implemented as a generator.

    :param f: function of one scalar to sample, should be positive,
              superflous kwargs are passed to it
    :param interval: the interval to sample from
    :param upper_bound: an upper bound to the function, optional
    :param seed: the seed for the rng, if not specified, the system
        time is used
    :param chunk_size: the size of the chunks of random numbers
        allocated per unit interval
    :yields: random nubers following the distribution of f
    :rtype: float
    """

    interval = _process_interval(interval)
    interval_length = (interval[1] - interval[0])

    if not upper_bound:
        upper_bound = find_upper_bound(f, interval, **kwargs)

    def allocate_random_chunk():
        return np.random.uniform([interval[0], 0], [interval[1], 1],
                            [chunk_size*interval_length, 2])

    while True:
        points = allocate_random_chunk()
        sample_points = points[:, 0] \
            [np.where(f(points[:, 0]) > points[:, 1]*upper_bound)]
        for point in sample_points:
            yield point
