# Author(s): Wojciech Reise
#
# Copyright (C) 2022 Inria
import numpy as np


def intersecting_lines(n_points, noise_std, center=0):
    """Sample two lines intersecting perpendicularly at center, in the unit cube.
    The sample is of size n_points in total and is corrupted by additive noise.
    :param n_points: number of points returned
    :type n_points: int
    :param noise_std: standard deviation of the noise
    :type noise_std: positive float
    :param center: point of intersection of the lines without noise
    :type center: np.ndarray of shape (1, 2) or 0.
    """
    n_ = n_points//2
    _n = n_points - n_

    t_, _t = np.linspace(0, 1, n_), np.linspace(0, 1, _n)
    upwards, downwards = np.stack([t_, t_], axis=1), np.stack([_t, -_t+1], axis=1)
    pts = np.concatenate([upwards, downwards], axis=0)
    pts += noise_std*np.random.randn(pts.shape[0], 2)

    pts += center - np.array([[0.5, 0.5]])
    return pts
