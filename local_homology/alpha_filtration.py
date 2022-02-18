# Author(s): Wojciech Reise
#
# Copyright (C) 2022 Inria

import numpy as np

from sklearn.metrics import pairwise_distances

from gudhi import RipsComplex, SimplexTree


## Coning Variant
def compute_local_homology_alpha(point_cloud, x0, epsilon, max_dimension, expand_boundary=True):
    """Compute the local homology of point_cloud at x0, with a pseudo alpha-filtration.

    """
    st = build_local_complex_alpha(point_cloud, x0, epsilon, max_dimension + 1,
                                   expand_boundary=expand_boundary)
    pd = get_persistence(st, max_dimension)
    return pd


def is_point_in_ball(point_cloud, point, epsilon, return_distances=True):
    """Determine which points from point_cloud are in the epsilon-ball centered at point.
    :param point_cloud: coordinates of points.
    :type point_cloud: np.ndarray of size (n_pts, dim)
    :param point: center of the ball that defines the neighborhood, does not have to be in point_cloud.
    :type point: np.ndarray of shape (1, dim)
    :param epsilon: size of the neighborhood to consider
    :type epsilon: float, positive
    :param return_distances: return distances if True
    :type return_distances: bool
    
    :returns: arrays specifying whether a given ball is, or not, in the complex.
    :
    """
    distance_from_point = pairwise_distances(point_cloud, point, metric="euclidean")
    distance_to_boundary = epsilon - distance_from_point
    is_in_ball = (distance_to_boundary >= 0)[:, 0]
    if return_distances:
        return is_in_ball, distance_to_boundary
    return is_in_ball


def build_local_complex_alpha(point_cloud, point, epsilon, max_dimension, expand_boundary=True):
    """Build a simplicial complex on vertices from
    point_cloud \cap B(point, \epsilon), along with an extra vertex
    representing the boundary.
    :param point_cloud: coordinates of points.
    :type point_cloud: np.ndarray of size (n_pts, dim)
    :param point: center of the ball that defines the neighborhood,
        does not have to be in point_cloud.
    :type point: np.ndarray of shape (1, dim)
    :param epsilon: size of the neighborhood to consider
    :type epsilon: float, positive
    :param max_dimension: max dimension to compute the homology of.
    :param expand_boundary: expand boundary inwards. Equivalent to including infinitely
        many points on the boundary.
    :type expand_boundary: bool
    """
    point_cloud_dist = pairwise_distances(point_cloud, metric="euclidean")
    is_in_ball, distance_to_boundary = is_point_in_ball(point_cloud, point,
                                                        epsilon, return_distances=True)
    
    point_cloud_dist_local = point_cloud_dist[is_in_ball][:, is_in_ball]
    local_boundary = distance_to_boundary[is_in_ball]

    if expand_boundary:
        local_boundary /= 2.

    point_cloud_local_stared = np.concatenate([point_cloud_dist_local, local_boundary.T],
                                              axis=0)
    point_cloud_local_stared = np.concatenate([point_cloud_local_stared,
                                               np.concatenate([local_boundary, np.array([[0]])],
                                                              axis=0)], axis=1)
    
    rc = RipsComplex(distance_matrix=point_cloud_local_stared)
    st = rc.create_simplex_tree(max_dimension=max_dimension)
    st.expansion(max_dimension)
    n_points_in_ball = np.sum(is_in_ball)
    star_vertex = n_points_in_ball  # indexing starts at 0!
    st.assign_filtration([star_vertex], -1)
    return st


def get_persistence(simplexTree, max_dimension):
    """Calculate the persistent homology of the abstract simplicial complex,
    filtering by positive values and dimensions.
    :param simplexTree: a simplcial complex, as returned by `build_local_complex`
    :type simplexTree: simplexTree
    :param max_dimension: max dimension of persistent homology to be returned.
    :type max_dimension: int.
    
    :returns: persistence diagram in the Gudhi format.
    :rtype: list of tuples (dim, (birth, death)).
    """
    pd = simplexTree.persistence()
    return [p for p in pd if (p[0] <= max_dimension) & (p[1][0] >= 0.)]
