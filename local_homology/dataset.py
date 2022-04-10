# Author(s): Wojciech Reise
#
# Copyright (C) 2022 Inria
import numpy as np
from scipy.spatial.transform import Rotation


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


def rectangle_with_hole(n_points, noise=0.0, return_rotation_matrix=False):
    """Sample points from a rectangle with a hole, embedded in 3d with a random rotation.
    The implementation is inspired by the Swiss hole from scikit-learn.
    :return X_in_3d: np.ndarray of shape (3, n_points)"""
    corners = np.array(
        [[np.pi * (1.5 + i), j * 7] for i in range(3) for j in range(3)]
    )
    corners = np.delete(corners, 4, axis=0)
    corner_index = np.random.choice(8, n_points)
    parameters = np.random.rand(2, n_points) * np.array([[np.pi], [7]])
    t, y = corners[corner_index].T + parameters

    X = np.stack([t, y, np.zeros(y.shape[0])], axis=1)

    random_rotation = Rotation.random(1).as_matrix()[0]
    X_in_3d = random_rotation @ (X.T)
    print(random_rotation.shape, X.shape, X_in_3d.shape)
    X_in_3d += noise * np.random.randn(3, n_points)
    if return_rotation_matrix:
        return X_in_3d.T, t, random_rotation
    return X_in_3d.T, t
