import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from monte_carlo import *


def plot_cubes(f, increments):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    cubes = generate_cubes(increments)
    for cube in cubes:
        x_range, y_range = cube
        xs, ys = np.mgrid[
            x_range[0] : x_range[1] : 0.01, y_range[0] : y_range[1] : 0.01,
        ]

        vol = get_integration_volume(cube)
        fun = f(xs, ys)

        ax.plot_surface(xs, ys, fun * vol * len(cubes), cmap=cm.coolwarm, linewidth=0)
        ax.plot_surface(xs, ys, fun, cmap=cm.coolwarm, linewidth=0)
        # ax.plot_surface(xs, ys, np.ones_like(fun) / vol, cmap=cm.coolwarm, linewidth=0)

    return fig, ax


def test_f(f, *args, **kwargs):
    res = integrate_vegas_nd(f, *args, **kwargs)
    fig, ax = plot_cubes(f, res.increment_borders)
    plt.show()
    return res
