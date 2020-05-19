"""
Simple monte carlo integration implementation.
Author: Valentin Boettcher <hiro@protagon.space>
"""

import os
import os.path
import numpy as np
import pathlib
import time
import numba
import collections
from typing import Union, List, Tuple, Callable, Iterator
from multiprocessing import Pool, cpu_count
import functools
from scipy.optimize import minimize_scalar, root, shgo, minimize
from dataclasses import dataclass
import utility
import joblib


def _process_interval(interval):
    interval = np.asarray(interval)
    if len(interval.shape) > 1:
        return np.array([_process_interval(i) for i in interval])

    assert len(interval) == 2, "An interval has two endpoints"

    a, b = interval
    if b < a:
        a, b = b, a

    return np.array([a, b])


@dataclass
class IntegrationResult:
    """
    A simple container class to hold the result, uncertainty and
    sampling size of the naive monte-carlo integration.
    """

    result: float
    sigma: float
    N: int

    @property
    def combined_result(self):
        """
        Get the result and accuracy combined as tuple.
        """

        return self.result, self.sigma


@utility.numpy_cache("cache")
def integrate(
    f,
    interval,
    epsilon=0.01,
    seed=None,
    num_points=None,
    adapt=True,
    return_sample=False,
    **kwargs,
) -> IntegrationResult:
    """Monte-Carlo integrates the function `f` over an interval.

    If the integrand is multi-dimensional it must accept an array of
    argument arrays. Think about your axes!

    :param f: function of one variable, kwargs are passed to it
    :param tuple interval: a 2-tuple of numbers, specifiying the
        integration range
    :param epsilon: desired accuracy
    :param seed: the seed for the rng, if not specified, the system
        time is used
    :param bool return_sample: wether to return the utilized sample
        of f as second return value

    :returns: the integration result

    :rtype: IntegrationResult

    """

    interval = _process_interval(interval)

    if len(interval.shape) == 1:
        return integrate(f, [interval], epsilon, seed, **kwargs)

    if seed:
        np.random.seed(seed)

    integration_volume = (interval[:, 1] - interval[:, 0]).prod()
    dimension = len(interval)

    # guess the correct N
    probe_points = np.random.uniform(
        interval[:, 0],
        interval[:, 1],
        (int(integration_volume * dimension * 100), dimension),
    ).T

    if num_points is None:
        prelim_std = integration_volume * f(*probe_points, **kwargs).std()
        epsilon = epsilon if prelim_std > epsilon else prelim_std / 10

        num_points = int((prelim_std / epsilon) ** 2 * 1.1 + 1)

    # now we iterate until we hit the desired epsilon
    while True:
        points = np.random.uniform(
            interval[:, 0], interval[:, 1], (num_points, len(interval))
        ).T

        sample = f(*points, **kwargs)
        integral = np.sum(sample) / num_points * integration_volume

        # the deviation gets multiplied by the square root of the interval
        # lenght, because it is the standard deviation of the integral.
        sample_std = np.std(sample) * integration_volume
        deviation = sample_std * np.sqrt(1 / (num_points - 1))

        if not adapt or deviation < epsilon:
            result = IntegrationResult(integral, deviation, num_points)

            return (result, sample) if return_sample else result

        # then we refine our guess, the factor 1.1
        num_points = int((sample_std / epsilon) ** 2 * 1.1)


def _negate(f):
    """A helper that multiplies the given function with -1."""

    @functools.wraps(f)
    def negated(*args, **kwargs):
        return -f(*args, **kwargs)

    return negated


def find_upper_bound(f, interval, **kwargs):
    """Find the upper bound of a function.

    :param f: function of one scalar and some kwargs that are passed
              on to it
    :param interval: interval to look in

    :returns: the upper bound of the function
    :rtype: float
    """

    upper_bound = minimize_scalar(
        lambda *args: -f(*args, **kwargs), bounds=interval, method="bounded"
    )
    if upper_bound.success:
        return -upper_bound.fun
    else:
        raise RuntimeError("Could not find an upper bound.")


def find_upper_bound_vector(f, interval, x0=None, tolerance=0.01):
    result = minimize(lambda x: -f(*x), x0=x0, bounds=interval, tol=tolerance)

    if not result.success:
        raise RuntimeError("Could not find an upper bound.", result)

    upper_bound = -result.fun
    return upper_bound


def sample_unweighted_vector(
    f, interval, seed=None, upper_bound=None, report_efficiency=False, total=1,
):
    dimension = len(interval)
    interval = _process_interval(interval)

    if seed:
        np.random.seed(seed)

    if upper_bound is None:
        upper_bound = find_upper_bound_vector(f, interval)

    def allocate_random_chunk():
        return np.random.uniform(
            [*interval[:, 0], 0], [*interval[:, 1], 1], [1, 1 + dimension],
        )

    total_points = 0
    total_accepted = 0

    while True:
        points = allocate_random_chunk()

        if report_efficiency:
            total_points += 1

        arg = points[:, 0:-1][0]
        if f(arg) >= points[:, -1] * upper_bound:
            if report_efficiency:
                total_accepted += 1

            yield (arg, total_accepted / total_points,) if report_efficiency else arg

    return


def sample_weighted_vector(f, interval, num_samples=1, seed=None):
    dimension = len(interval)
    interval = _process_interval(interval)

    if seed:
        np.random.seed(seed)

    sample_points = np.random.uniform(
        interval[:, 0], interval[:, 1], [num_samples, dimension],
    )

    return sample_points, np.array([f(sample) for sample in sample_points])


def sample_unweighted(
    f,
    interval,
    upper_bound=None,
    seed=None,
    chunk_size=100,
    report_efficiency=False,
    **kwargs,
):
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
    interval_length = interval[1] - interval[0]

    if seed:
        np.random.seed(seed)

    upper_bound_fn, upper_bound_integral, upper_bound_integral_inverse = (
        None,
        None,
        None,
    )
    # i know....

    if not upper_bound:
        upper_bound_value = find_upper_bound(f, interval, **kwargs)

        def upper_bound_fn(x):
            return upper_bound_value

        def upper_bound_integral(x):
            return upper_bound_value * x

        def upper_bound_integral_inverse(y):
            return y / upper_bound_value

    elif len(upper_bound) == 2:
        upper_bound_fn, upper_bound_integral = upper_bound

        def upper_inv(points):  # not for performance right now...
            return np.array(
                [
                    root(
                        lambda y: upper_bound_integral(y) - x, x0=0, jac=upper_bound_fn
                    ).x
                    for x in points
                ]
            ).T

        upper_bound_integral_inverse = upper_inv

    elif len(upper_bound) == 3:
        upper_bound_fn, upper_bound_integral, upper_bound_integral_inverse = upper_bound
    else:
        raise ValueError("The upper bound must be `None` or a three element sequence!")

    def allocate_random_chunk():
        return np.random.uniform(
            [upper_bound_integral(interval[0]), 0],
            [upper_bound_integral(interval[1]), 1],
            [int(chunk_size * interval_length), 2],
        )

    total_points = 0
    total_accepted = 0

    while True:
        points = allocate_random_chunk()
        points[:, 0] = upper_bound_integral_inverse(points[:, 0])
        sample_points = points[:, 0][
            np.where(f(points[:, 0]) > points[:, 1] * upper_bound_fn(points[:, 0]))
        ]

        if report_efficiency:
            total_points += points.size
            total_accepted += sample_points.size

        for point in sample_points:
            yield (point, total_accepted / total_points) if report_efficiency else point


@dataclass
class VegasIntegrationResult(IntegrationResult):
    """
    A simple container class to hold the result, uncertainty and
    sampling size of the naive monte-carlo integration.

    The ascertained increment borders, as well as the total sample
    size are available as well.
    """

    increment_borders: np.ndarray
    vegas_iterations: int


def reshuffle_increments(
    integral_steps,
    integral,
    interval_lengths,
    num_increments,
    increment_borders,
    alpha,
    K,
):

    # alpha controls the convergence speed
    μ = np.abs(integral_steps) / integral
    new_increments = (K * ((μ - 1) / (np.log(μ))) ** alpha).astype(int)
    group_size = new_increments.sum() / num_increments
    new_increment_borders = np.empty_like(increment_borders)

    # this whole code does a very simple thing: it eats up
    # sub-increments until it has `group_size` of them
    i = 0  # position in increment count list
    j = 0  # position in new_incerement_borders

    # the number of sub-increments still available
    rest = new_increments[0]

    # the number of sub-increments needed to fill one increment
    head = group_size

    # the current position in the interval relative to its
    # beginning
    current = 0

    while i < num_increments and (j < (num_increments - 1)):
        if new_increments[i] == 0:
            i += 1
            rest = new_increments[i]

        current_increment_size = interval_lengths[i] / new_increments[i]

        if head <= rest:
            current += head * current_increment_size
            new_increment_borders[j] = current
            rest -= head
            head = group_size
            j += 1

        else:
            current += rest * current_increment_size
            head -= rest
            i += 1
            rest = new_increments[i]

    return new_increment_borders


def integrate_vegas(
    f,
    interval,
    seed=None,
    num_increments=5,
    epsilon=1e-3,
    increment_epsilon=1e-2,
    alpha=1.5,
    acumulate=True,
    vegas_point_density=1000,
    **kwargs,
) -> VegasIntegrationResult:
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
    :param increment_epsilon: the breaking condition, if the magnitude
        of the difference between the increment positions in
        subsequent iterations does not change more then epsilon the
        computation is considered to have converged
    :param alpha: controls the the convergence speed, should be
        between 1 and 2 (the lower the faster)
    :returns: the intregal, the standard deviation, an array of
              increment borders which can be used in subsequent
              sampling
    :rtype: tuple
    """

    interval = _process_interval(interval)
    interval_length = interval[1] - interval[0]

    if seed:
        np.random.seed(seed)

    # no clever logic is being used to define the vegas iteration
    # sample density for the sake of simplicity
    points_per_increment = int(vegas_point_density * interval_length / num_increments)

    # start with equally sized intervals
    interval_borders = np.linspace(*interval, num_increments + 1, endpoint=True)

    def evaluate_integrand(interval_borders, interval_lengths, samples_per_increment):
        intervals = np.array((interval_borders[:-1], interval_borders[1:]))
        sample_points = np.random.uniform(
            *intervals, (samples_per_increment, num_increments)
        ).T

        weighted_f_values = f(sample_points, **kwargs) * interval_lengths[:, None]

        # the mean here has absorbed the num_increments
        integral_steps = weighted_f_values.mean(axis=1)
        integral = integral_steps.sum()
        variance = (
            (f(sample_points, **kwargs).std(axis=1) * interval_lengths) ** 2
        ).sum() / (samples_per_increment - 1)
        return integral, integral_steps, variance

    K = num_increments * 1000
    increment_borders = interval_borders[1:-1] - interval_borders[0]

    integrals = []
    variances = []

    vegas_iterations, integral, variance = 0, 0, 0
    while True:
        vegas_iterations += 1
        interval_lengths = interval_borders[1:] - interval_borders[:-1]
        integral, integral_steps, variance = evaluate_integrand(
            interval_borders, interval_lengths, points_per_increment
        )

        integrals.append(integral)
        variances.append(variance)

        # it is debatable to pass so much that could be recomputed...
        new_increment_borders = reshuffle_increments(
            integral_steps,
            integral,
            interval_lengths,
            num_increments,
            increment_borders,
            alpha,
            K,
        )

        interval_borders[1:-1] = interval_borders[0] + increment_borders
        if (
            np.linalg.norm(increment_borders - new_increment_borders)
            < increment_epsilon
        ):
            break

        increment_borders = new_increment_borders

    interval_lengths = interval_borders[1:] - interval_borders[:-1]

    # brute force increase of the sample size
    if np.sqrt(variance) >= epsilon:
        tick = 3
        while True:
            integral, _, variance = evaluate_integrand(
                interval_borders, interval_lengths, points_per_increment
            )

            integrals.append(integral)
            variances.append(variance)

            if np.sqrt(variance) <= epsilon:
                break

            # adaptive scaling of sample size incrementation
            points_per_increment += int(
                1000 ** np.log(tick) * interval_length / num_increments
            )

            tick += 2

    # as a bonus, we utilize all prior integration results
    if acumulate:
        integrals = np.array(integrals)
        variances = np.array(variances)
        integral = np.sum(integrals ** 3 / variances ** 2) / np.sum(
            integrals ** 2 / variances ** 2
        )
        variance = 1 / np.sqrt(np.sum(integrals ** 2 / variances ** 2)) * integral

    return VegasIntegrationResult(
        integral,
        np.sqrt(variance),
        points_per_increment * num_increments,
        interval_borders,
        vegas_iterations,
    )


def sample_stratified(
    f, increment_borders, seed=None, chunk_size=100, report_efficiency=False, **kwargs
):
    """Samples a distribution proportional to f by hit and miss.
    Implemented as a generator.

    :param f: function of one scalar to sample, should be positive,
              superflous kwargs are passed to it
    :param interval: the interval to sample from
    :param seed: the seed for the rng, if not specified, the system
        time is used
    :param chunk_size: the size of the chunks of random numbers
        allocated per unit interval
    :yields: random nubers following the distribution of f
    """

    increment_count = increment_borders.size - 1
    increment_lenghts = increment_borders[1:] - increment_borders[:-1]
    weights = increment_count * increment_lenghts
    increment_chunk = int(chunk_size / increment_count)
    chunk_size = increment_chunk * increment_count

    upper_bound = np.array(
        [
            find_upper_bound(
                lambda x: f(x, **kwargs) * weight, [left_border, right_border]
            )
            for weight, left_border, right_border in zip(
                weights, increment_borders[:-1], increment_borders[1:]
            )
        ]
    ).max()

    total_samples = 0
    total_accepted = 0

    while True:
        increment_samples = np.random.uniform(
            increment_borders[:-1],
            increment_borders[1:],
            [increment_chunk, increment_count],
        )
        increment_y_samples = np.random.uniform(
            0, 1, [increment_chunk, increment_count]
        )
        f_weighted = f(increment_samples) * weights  # numpy magic at work here
        mask = f_weighted > increment_y_samples * upper_bound

        if report_efficiency:
            total_samples += chunk_size
            total_accepted += np.count_nonzero(mask)

        for point in increment_samples[mask]:
            yield (
                point,
                total_accepted / total_samples,
            ) if report_efficiency else point


class SamplingWorker:
    def __init__(self, num_samples, **kwargs):
        self.num_samples = num_samples
        self.kwargs = kwargs

    def draw_sample(self):
        global _FUN
        seed = int.from_bytes(os.urandom(2), "big")
        return sample_unweighted_array(self.num_samples, _FUN, **self.kwargs, seed=seed)


_FUN = None


@utility.numpy_cache("cache")
def sample_unweighted_array(
    num,
    f,
    interval=None,
    increment_borders=None,
    cubes=None,
    report_efficiency=False,
    proc=None,
    status_path=None,
    **kwargs,
):
    """Sample `num` elements from a distribution.  The rest of the
    arguments is analogous to `sample_unweighted`.
    """
    global _FUN

    sample_arr = None
    eff = None

    # some multiprocessing magic
    if proc is not None:
        if isinstance(proc, str) and proc == "auto":
            proc = cpu_count()

        result = None
        num_per_worker = int(num / proc + 1)
        _FUN = f  # there is no other way :(

        workers = [
            SamplingWorker(
                num_samples=num_per_worker,
                interval=interval,
                increment_borders=increment_borders,
                report_efficiency=report_efficiency,
                status_path=status_path,
                proc=None,
                cache=None,
                cubes=cubes,
                **kwargs,
            )
            for _ in range(proc)
        ]

        with Pool(proc) as p:
            result = p.map(SamplingWorker.draw_sample, workers)
            result = np.array(result)

        if report_efficiency:
            eff = result[:, 1].mean()
            sample_arr = np.concatenate(result[:, 0])[0:num]
        else:
            sample_arr = np.concatenate(result)[0:num]

    else:
        samples = None

        if interval is not None:
            interval = np.array(interval)
            vectorized = len(interval.shape) > 1
            sample_arr = np.empty((num, interval.shape[0]) if vectorized else num)
            if len(interval.shape) > 1:
                samples = sample_unweighted_vector(
                    f,
                    interval,
                    report_efficiency=report_efficiency,
                    total=num,
                    **kwargs,
                )
            else:
                if "chunk_size" not in kwargs:
                    kwargs["chunk_size"] = num * 10

                samples = sample_unweighted(
                    f, interval, report_efficiency=report_efficiency, **kwargs
                )
        elif cubes is not None:
            sample_arr = np.empty((num, len(cubes[0][0])))
            samples = sample_stratified_vector(
                f, cubes, report_efficiency=report_efficiency, **kwargs,
            )
        elif increment_borders is not None:
            sample_arr = np.empty(num)
            samples = sample_stratified(
                f,
                increment_borders=increment_borders,
                report_efficiency=report_efficiency,
                **kwargs,
            )
        else:
            raise TypeError("Neiter interval nor increment_borders or cubes specified!")

        fifo = None
        if status_path:
            if not pathlib.Path(status_path).exists():
                os.mkfifo(status_path)
            fifo = open(status_path, "w")

        start_time = time.monotonic()
        last_time = start_time

        for i, sample in zip(range(num), samples):
            if fifo:
                this_time = time.monotonic()
                if this_time - last_time > 1:
                    δ = this_time - start_time
                    total = i + 1
                    speed = total / δ

                    try:
                        fifo.write(
                            f"{os.getpid()}: {i:<10} "
                            f"{num-total:<10} {speed:4.1f} 1/sec "
                            f"{(total / num) * 100:2.1f}%"
                            f"{(num-total)/speed/60:4.0f} min\n"
                        )
                        fifo.flush()
                    except BrokenPipeError:
                        pass

                    last_time = this_time  # after the work is before the work is ...

            if report_efficiency:
                sample_arr[i], _ = sample
            else:
                sample_arr[i] = sample

        eff = next(samples)[1] if report_efficiency else None

        fifo and fifo.close()

    return (sample_arr, eff) if report_efficiency else sample_arr


# Geezus! Such code dublication..., but I don't want to break
# backwards compatibility.

CubeType = List[Tuple[np.ndarray, float, float]]


@dataclass
class VegasIntegrationResultNd(VegasIntegrationResult):
    """A simple container class to hold the result, uncertainty and
    sampling size of the multi-dimensional vegas integration.

    The ascertained increment borders, the hypercubes with their
    respective maxima as well as the total sample size are available
    as well.
    """

    increment_borders: np.ndarray
    cubes: CubeType
    vegas_iterations: int


@utility.numpy_cache("cache")
def integrate_vegas_nd(
    f,
    intervals,
    seed=None,
    num_increments=5,
    epsilon=1e-3,
    increment_epsilon=1e-2,
    alpha=1.5,
    vegas_point_density=1000,
    proc="auto",
    num_points_per_cube=None,
    vegas_step_f=None,
    **kwargs,
) -> VegasIntegrationResult:
    """Integrate the given function (in n-dimensions) with the vegas
    algorithm to reduce variance.  This implementation follows the
    description given in JOURNAL OF COMPUTATIONAL 27, 192-203 (1978).

    :param f: function of one variable, kwargs are passed to it
    :param tuple interval: a 2-tuple of numbers, specifiying the
        integration range
    :param seed: the seed for the rng, if not specified, the system
        time is used
    :param num_increments: the number increments in which to divide
        the interval
    :param point_density: the number of random points per unit
        interval
    :param increment_epsilon: the breaking condition, if the magnitude
        of the difference between the increment positions in
        subsequent iterations does not change more then epsilon the
        computation is considered to have converged
    :param alpha: controls the the convergence speed, should be
        between 1 and 2 (the lower the faster)

    :returns: the intregal, the standard deviation, an array of
              increment borders which can be used in subsequent
              sampling

    :rtype: tuple

    """
    intervals = np.asarray(_process_interval(intervals))
    ndim = len(intervals)
    integration_volume = (intervals[:, 1] - intervals[:, 0]).prod()

    if not isinstance(num_increments, collections.Iterable):
        num_increments = np.ones(ndim).astype(int) * num_increments
    else:
        num_increments = np.asarray(num_increments)

    num_cubes = num_increments.prod()

    if seed:
        np.random.seed(seed)

    if proc == "auto":
        proc = cpu_count()

    # no clever logic is being used to define the vegas iteration
    # sample density for the sake of simplicity
    points_per_cube = num_points_per_cube or int(
        vegas_point_density * integration_volume / num_cubes
    )

    # start with equally sized intervals
    increment_borders = [
        np.linspace(*interval, num_increments[i] + 1, endpoint=True)
        for i, interval in enumerate(intervals)
    ]

    vegas_step_f = vegas_step_f or f

    def evaluate_stripe(interval_borders):
        cubes = generate_cubes(interval_borders)

        result = 0
        for cube in cubes:
            vol = get_integration_volume(cube)
            points = np.random.uniform(
                cube[:, 0], cube[:, 1], (points_per_cube, ndim)
            ).T
            result += ((vegas_step_f(*points, **kwargs) * vol) ** 2).sum()

        return np.sqrt(result)

    def generate_integral_steps(interval_borders, dimension):
        borders = interval_borders[dimension]
        stripes = np.array([borders[:-1], borders[1:]]).T
        ms = np.empty(len(stripes))

        for i, stripe in enumerate(stripes):
            interval_borders[dimension] = stripe

            ms[i] = evaluate_stripe(interval_borders)

        # * num_increments gets normalized away
        return ms / ms.sum()

    K = 1000 * num_increments

    vegas_iterations, integral, variance = 0, 0, 0

    while True:
        vegas_iterations += 1
        new_increment_borders = []
        remainder = 0

        for dim in range(ndim):
            increment_weights = generate_integral_steps(increment_borders.copy(), dim)
            nonzero_increments = increment_weights[increment_weights > 0]
            if nonzero_increments.size > 0:
                remainder += nonzero_increments.max() - nonzero_increments.min()

            new_borders = reshuffle_increments_nd(
                increment_weights,
                num_increments[dim],
                increment_borders[dim],
                alpha,
                K[dim],
            )
            new_increment_borders.append(new_borders)

        remainder /= ndim
        print(remainder)
        increment_borders = new_increment_borders
        if abs(remainder) < increment_epsilon:
            break

    # brute force increase of the sample size
    cubes = generate_cubes(increment_borders)
    volumes = [get_integration_volume(cube) for cube in cubes]
    cube_samples = [np.empty(0) for _ in cubes]
    cube_sample_points = [np.empty((0, ndim)) for _ in cubes]
    total_points_per_cube = 0
    maxima = np.empty((num_cubes, 1 + ndim))
    integrals = np.empty(num_cubes)

    @joblib.delayed
    def evaluate_integrand(cube):
        points = np.random.uniform(cube[:, 0], cube[:, 1], (points_per_cube, ndim)).T
        return f(*points, **kwargs), points

    with joblib.Parallel(n_jobs=proc) as parallel:
        while True:
            integral = variance = 0
            total_points_per_cube += points_per_cube
            samples = parallel([evaluate_integrand(cube) for cube in cubes])

            for (sample, points), vol, i in zip(samples, volumes, range(num_cubes)):
                # let's re-use the samples from earlier runs
                cube_samples[i] = np.concatenate([cube_samples[i], sample]).T
                cube_sample_points[i] = np.concatenate(
                    [cube_sample_points[i], points.T]
                )
                points = cube_samples[i] * vol
                curr_integral = points.mean()
                integral += curr_integral
                variance += (
                    cube_samples[i].var() * (vol ** 2) / (total_points_per_cube - 1)
                )

                max_index = cube_samples[i].argmax()
                maxima[i] = [
                    cube_samples[i][max_index],
                    *cube_sample_points[i][max_index],
                ]

                integrals[i] = curr_integral

            if np.sqrt(variance) <= epsilon:
                break

            # adaptive scaling of sample size incrementation
            points_per_cube *= int((variance / epsilon ** 2) * 1.5)
            print(points_per_cube, integral, variance)

    return VegasIntegrationResultNd(
        integral,
        np.sqrt(variance),
        total_points_per_cube * num_cubes,
        increment_borders,
        vegas_iterations,
        [
            (cube, maximum, integral)
            for cube, maximum, integral in zip(cubes, maxima, integrals)
        ],
    )


def generate_cubes(interval_borders):
    """Given an array of interval borders, return a list of hypercube
    edges to fill the whole volume."""
    intervals = [np.array([border[:-1], border[1:]]).T for border in interval_borders]
    ndim = len(intervals)
    axis = np.arange(ndim)

    mesh_axes = [np.arange(len(ax_intervals)) for ax_intervals in intervals]
    grid = np.array(np.meshgrid(*mesh_axes)).T
    grid = grid.reshape(int(grid.size / ndim), ndim).astype(int)

    return np.array(
        [[intervals[n][i] for (n, i) in zip(axis, indices)] for indices in grid]
    )


def get_integration_volume(interval):
    return (interval[:, 1] - interval[:, 0]).prod()


def reshuffle_increments_nd(
    μ, num_increments, increment_borders, alpha, K,
):
    if num_increments == 1:
        return increment_borders

    # alpha controls the convergence speed
    new_increments = np.empty_like(μ)

    if np.any(μ == 1):
        return increment_borders

    mask = μ > 0
    μ = μ[mask]
    new_increments[mask] = (K * ((μ - 1) / (np.log(μ))) ** alpha).astype(int)
    new_increments[np.logical_not(mask)] = 1

    group_size = int(new_increments.sum() / num_increments)

    new_increment_borders = np.empty_like(increment_borders[1:-1])
    interval_lengths = increment_borders[1:] - increment_borders[:-1]

    # this whole code does a very simple thing: it eats up
    # sub-increments until it has `group_size` of them
    i = 0  # position in increment count list
    j = 0  # position in new_incerement_borders

    # the number of sub-increments still available
    rest = new_increments[0]

    # the number of sub-increments needed to fill one increment
    head = group_size

    # the current position in the interval relative to its
    # beginning
    current = 0
    while i < num_increments and (j < (num_increments - 1)):
        if new_increments[i] == 0:
            i += 1
            rest = new_increments[i]

        current_increment_size = interval_lengths[i] / new_increments[i]

        if head <= rest:
            current += head * current_increment_size
            new_increment_borders[j] = current
            rest -= head
            head = group_size
            j += 1

        else:
            current += rest * current_increment_size
            head -= rest
            i += 1
            rest = new_increments[i]

    return (
        np.array(
            [0, *new_increment_borders, increment_borders[-1] - increment_borders[0]]
        )
        + increment_borders[0]
    )


def sample_stratified_vector(
    f: Callable[..., Union[float, int]],
    cubes: CubeType,
    chunk_size: int = 1,
    seed: float = None,
    report_efficiency: bool = False,
    overestimate_factor: float = 1,
) -> Iterator[Union[np.ndarray, Tuple[np.ndarray, float]]]:
    """Sample a distribution `f` using stratified sampling and hit-or-miss
    in the subvolumes `cubes` which contain information about the
    total probability per cube and the upper-bound per cube.

    :param f: the distribution, function of multiple parameters returning a number
    :param cubes: a list of subvolumes,
      elements of the form [cube borders, maximum, integral]
    :param chunk_size: the number of elements to take consecutively per subvolume,
      should be 1 for most purposes, the higher this parameter the faster/unreliable
      the sampling gets
    :param seed: the seed of the RNG, if none
    :param report_efficiency: Wether to yield the sampling efficiency
    :param overestimate_factor: by how much to over_estimeate the maximum
    :yields: an ndarray containing a sample and optionally a the sampling efficiency
    """

    ndim = len(cubes[0][0])

    if seed:
        np.random.seed(seed)

    total_points = 0
    total_accepted = 0

    cubes = [cube for cube in cubes if cube[2] > 0]  # filter out cubes with zero weight
    weights = np.array([weight for _, _, weight in cubes])
    integral = weights.sum()
    weights = np.cumsum(weights / integral)

    maxima = np.array(
        [
            find_upper_bound_vector(
                f, cube[0], cube[1][1:], tolerance=overestimate_factor - 1
            )
            for cube in cubes
        ]
    )

    # compiling this function results in quite a speed-up
    @numba.jit(numba.int32(), nopython=True)
    def find_next_index():
        cube_index = 0
        rand = np.random.uniform(0, 1)
        for weight in weights:
            if weight > rand:
                break

            cube_index += 1

        return cube_index

    while True:
        cube_index = find_next_index()

        counter = 0
        cube, _, _ = cubes[cube_index]
        maximum = maxima[cube_index]

        while counter < chunk_size:
            points = np.random.uniform(
                [*cube[:, 0], 0], [*cube[:, 1], 1], [1, 1 + ndim],
            )

            if report_efficiency:
                total_points += 1

            args = points[0][:-1]
            if f(*args) >= points[0][-1] * maximum * overestimate_factor:
                if report_efficiency:
                    total_accepted += 1

                yield (
                    args.T,
                    total_accepted / total_points,
                ) if report_efficiency else args.T

                counter += 1

    return


def estimate_stratified_efficiency(cubes: CubeType) -> float:
    """Estimate the sampling efficiency of stratified sampling based on
the data about subvolumes.

    :param cubes: the hypercubes used in the stratified sampling
    :returns: the estimated sampling efficiency in percent
    """

    weights = np.array([weight for _, _, weight in cubes if weight > 0])
    maxima = np.array([maximum[0] for _, maximum, weight in cubes if weight > 0])
    volumes = np.array(
        [get_integration_volume(cube) for cube, _, weight in cubes if weight > 0]
    )

    return np.sum(weights) / np.sum((volumes * maxima)) * 100
