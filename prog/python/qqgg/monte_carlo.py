"""
Simple monte carlo integration implementation.
Author: Valentin Boettcher <hiro@protagon.space>
"""
import numpy as np
from scipy.optimize import minimize_scalar, root

def _process_interval(interval):
    assert len(interval) == 2, 'An interval has two endpoints'

    a, b = interval
    if b < a:
        a, b = b, a

    return [a, b]


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

    # the deviation gets multiplied by the square root of the interval
    # lenght, because it is the standard deviation of the integral.
    deviation = np.std(sample)*np.sqrt(1/(num_points - 1))*interval_length

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
                    chunk_size=100, report_efficiency=False, **kwargs):
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


    upper_bound_fn, upper_bound_integral, upper_bound_integral_inverse = None, None, None
    # i know....

    if not upper_bound:
        upper_bound_value = find_upper_bound(f, interval, **kwargs)
        def upper_bound_fn(x): return upper_bound_value
        def upper_bound_integral(x): return upper_bound_value*x
        def upper_bound_integral_inverse(y): return y/upper_bound_value

    elif len(upper_bound) == 2:
        upper_bound_fn, upper_bound_integral =\
            upper_bound

        def upper_inv(points):  # not for performance right now...
            return np.array([root(lambda y: upper_bound_integral(y) - x, x0=0,
                            jac=upper_bound_fn).x for x in points]).T

        upper_bound_integral_inverse = upper_inv

    elif len(upper_bound) == 3:
        upper_bound_fn, upper_bound_integral, upper_bound_integral_inverse =\
            upper_bound
    else:
        raise ValueError('The upper bound must be `None` or a three element sequence!')

    def allocate_random_chunk():
        return np.random.uniform([upper_bound_integral(interval[0]), 0],
                            [upper_bound_integral(interval[1]), 1],
                            [int(chunk_size*interval_length), 2])

    total_points = 0
    total_accepted = 0

    while True:
        points = allocate_random_chunk()
        points[:, 0] = upper_bound_integral_inverse(points[:, 0])
        sample_points = points[:, 0] \
            [np.where(f(points[:, 0]) > \
                      points[:, 1]*upper_bound_fn(points[:, 0]))]

        if report_efficiency:
            total_points += points.size
            total_accepted += sample_points.size

        for point in sample_points:
            yield (point, total_accepted/total_points) \
                if report_efficiency else point

def sample_unweighted_array(num, *args, report_efficiency=False, **kwargs):
    """Sample `num` elements from a distribution.  The rest of the
    arguments is analogous to `sample_unweighted`.
    """

    sample_arr = np.empty(num)
    samples = sample_unweighted(*args, report_efficiency=report_efficiency,
                                **kwargs)

    for i, sample in zip(range(num), samples):
        if report_efficiency:
            sample_arr[i], _ = sample
        else:
            sample_arr[i] = sample

    return (sample_arr, next(samples)[1]) if report_efficiency else sample_arr


def integrate_vegas(f, interval, seed=None, num_increments=5,
                  point_density=1000, epsilon=1e-2, alpha=1.5, **kwargs):
    """Integrate the given function (in one dimension) with the vegas
    algorithm to reduce variance.  This implementation follows the
    description given in JOURNAL OF COMPUTATIONAL 27, 192-203 (1978).

    All iterations contribute to the final result.

    :param f: function of one variable, kwargs are passed to it
    :param tuple interval: a 2-tuple of numbers, specifiying the
        integration range
    :param seed: the seed for the rng, if not specified, the system
        time is used
    :param num_increments: the number increments in which to divide
        the interval
    :param point_density: the number of random points per unit
        interval
    :param epsilon: the breaking condition, if the magnitude of the
        difference between the increment positions in subsequent
        iterations does not change more then epsilon the computation
        is considered to have converged
    :param alpha: controls the the convergence speed, should be
        between 1 and 2

    :returns: the intregal, the standard deviation, an array of
              increment borders which can be used in subsequent
              sampling

    :rtype: tuple
    """

    Interval = _process_interval(interval)
    interval_length = (interval[1] - interval[0])

    if seed:
        np.random.seed(seed)

    # start with equally sized intervals
    interval_borders = np.linspace(*interval, num_increments + 1, endpoint=True)

    points_per_increment = int(point_density*interval_length/num_increments)
    total_points = points_per_increment*num_increments

    def evaluate_integrand(interval_borders, interval_lengths):
        intervals = np.array((interval_borders[:-1], interval_borders[1:]))
        sample_points = np.random.uniform(*intervals,
                                          (points_per_increment, num_increments)).T

        weighted_f_values = f(sample_points, **kwargs)*interval_lengths[:, None]

        # here the num_increments don't cancel
        weighted_f_squared_values = (f(sample_points, **kwargs) \
            *interval_lengths[:, None])**2*num_increments

        integral_steps = weighted_f_values.mean(axis=1)
        integral = integral_steps.sum()
        variance = 1/(total_points - 1)\
            *(weighted_f_squared_values.mean(axis=1).sum() - integral**2)

        return integral_steps.sum(), integral_steps, variance

    K = num_increments*1000
    increment_borders = interval_borders[1:-1] - interval_borders[0]

    integrals = []
    variances = []

    while True:
        interval_lengths = interval_borders[1:] - interval_borders[:-1]
        integral, integral_steps, variance = \
            evaluate_integrand(interval_borders, interval_lengths)

        integrals.append(integral)
        variances.append(variance)

        # alpha controls the convergence speed
        μ = np.abs(integral_steps)/integral
        new_increments = (K*((μ - 1)/(np.log(μ)))**alpha).astype(int)
        group_size = new_increments.sum()/num_increments

        new_increment_borders = np.empty_like(increment_borders)

        # this whole code does a very simple thing: it eats up
        # sub-increments until it has `group_size` of them
        i = 0  # position in increment count list
        j = 0  # position in new_incerement_borders
        rest = new_increments[0]  # the number of sub-increments still available
        head = group_size  # the number of sub-increments needed to
                           # fill one increment
        current = 0  # the current position in the interval relative
                     # to its beginning

        while i < num_increments and (j < (num_increments - 1)):
            if new_increments[i] == 0:
                i += 1
                rest = new_increments[i]

            current_increment_size = interval_lengths[i]/new_increments[i]

            if head <= rest:
                current += head*current_increment_size
                new_increment_borders[j] = current
                rest -= head
                head = group_size
                j += 1

            else:
                current += rest*current_increment_size
                head -= rest
                i += 1
                rest = new_increments[i]

        interval_borders[1:-1] = interval_borders[0] + increment_borders
        if np.linalg.norm(increment_borders - new_increment_borders) < epsilon:

            integrals = np.array(integrals)
            variances = np.array(variances)
            integral = np.sum(integrals**3/variances**2) \
                /np.sum(integrals**2/variances**2)
            mean_variance = 1/np.sqrt(np.sum(integrals**2/variances**2))*integral
            return integral, np.sqrt(mean_variance), interval_borders

        increment_borders = new_increment_borders
