# Author(s): Wojciech Reise
#
# Copyright (C) 2022 Inria
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle


def plot_point_cloud(point_cloud):
    """Plot points in the plane."""
    plt.scatter(point_cloud[:, 0], point_cloud[:, 1])


def plot_disc(center, radius, **kwargs):
    """Plot a disc of radius `radius` centered at `center`."""
    alpha = kwargs.pop("alpha", 0.2)
    ax = plt.gca()
    ax.add_patch(plt.Circle(center, radius, alpha=alpha, **kwargs))


def plot_rectangle(v, **kwargs):
    alpha = kwargs.pop("alpha", 0.2)
    ax = plt.gca()
    (_, y_max), (x_min, _) = ax.get_ylim(), ax.get_xlim()
    ax.add_patch(Rectangle((x_min, v), width=v-x_min, height=y_max-v, alpha=alpha, **kwargs))