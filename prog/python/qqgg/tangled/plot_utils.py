import matplotlib.pyplot as plt

def draw_histo(points, xlabel, bins=20):
    heights, edges = np.histogram(points, bins)
    centers = (edges[1:] + edges[:-1])/2
    deviations = np.sqrt(heights)

    fig, ax = set_up_plot()
    ax.errorbar(centers, heights, deviations, linestyle='none', color='orange')
    ax.step(edges,  [heights[0], *heights], color='#1f77b4')

    ax.set_xlabel(xlabel)
    ax.set_xlim([points.min(), points.max()])
    return fig, ax
