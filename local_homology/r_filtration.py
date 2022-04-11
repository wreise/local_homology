# Author(s): Wojciech Reise
#
# Copyright (C) 2022 Inria

from typing import List, Tuple

import gudhi.simplex_tree
import numpy as np
from gudhi import SimplexTree
from matplotlib import pyplot as plt


def compute_local_homology_r(point_cloud: np.ndarray, x0: np.ndarray, alpha: float,
                             max_dimension: int, max_r: float = None) -> List[Tuple]:
    """Calculate the local persistent homology of `point_cloud` at point `x_0`,
    with scale parameter alpha. Computed using the approximation technique
    from [1]. **Note:** the minimal dimension with non-trivial
    features is 1.
    :param point_cloud: the point-cloud to compute the homology of.
    :type point_cloud:  ndarray of shape (n_points, d).
    :param x0: the point at which to compute the local homology.
    :type x0: ndarray of shape (1, d).
    :param alpha: construct a Rips with 2*alpha.
    :type alpha: positive float.
    :param max_dimension: max dimension of homology that you want to see.
    :type max_dimension: int

    :returns: list of tuples (dim, (b,d)) of the persistence diagram.

    1. Skraba, P. & Wang, B. Approximating Local Homology from Samples.
       In Proceedings of the Twenty-Fifth Annual ACM-SIAM Symposium on Discrete
       Algorithms 174–192 (Society for Industrial and Applied Mathematics, 2014).
       doi:10.1137/1.9781611973402.13.
    """
    St = get_filtered_complex(point_cloud, x0, alpha, max_dimension, max_r)
    # For H_k, by corollary, it is H_{k-1}, so you need to build up to dim=k.
    St.extend_filtration()

    dgm_ordinary_minus, _, _, _ = St.extended_persistence()
    dgm_relative = apply_duality_ordinary_to_relative(dgm_ordinary_minus)

    return dgm_relative


def get_filtered_complex(point_cloud: np.ndarray, x0: np.ndarray, alpha: float,
                         max_dimension: int, max_r: float = None) -> gudhi.simplex_tree.SimplexTree:
    """Construct the VietorisRips complex from `point_cloud` with scale alpha.
    Filter by the function $x \mapsto - d(x, x_0)$ and expand until
    dimension max_dimension.

    :param point_cloud: the point-cloud to build the complex of.
    :type point_cloud: ndarray of shape (n_points, d).
    :param x0: the point at which to compute the local homology.
    :type x0: ndarray of shape (1, d).
    :param alpha: construct a Rips with 2*alpha.
    :type alpha: positive float.
    :param max_dimension: maximal simplex dimension.
    :type alpha: int.
    :param max_r: exclude points which are too far from x_0
    :type max_r: float.

    :returns: gudhi.SimplexTree
    """
    dist_to_x0 = np.linalg.norm(point_cloud - x0, ord=2, axis=1)
    vertex_value = -dist_to_x0
    if max_r is not None:
        is_not_too_far = vertex_value >= - max_r
        vertex_value, point_cloud = vertex_value[is_not_too_far], point_cloud[is_not_too_far]
    pairwise_distances = np.linalg.norm(point_cloud[None, :, :] - point_cloud[:, None, :],
                                        ord=2, axis=-1)

    is_edge_in_complex = pairwise_distances < alpha

    St = SimplexTree()
    for i1, val in enumerate(vertex_value):
        St.insert([i1], val)

    for (i1, i2), is_in_complex in np.ndenumerate(is_edge_in_complex):
        if is_in_complex:
            edge_value = max(vertex_value[i1], vertex_value[i2])
            St.insert([i1, i2], edge_value)

    St.expansion(max_dimension)
    return St


def apply_duality_ordinary_to_relative(dgm_ordinary):
    """Apply the Symmetry Corollary, transforming the ordinary diagram of $-f$
    to the relative diagram of $f$.

    :param dgm_ordinary: represents a persistence diagram.
    :type dgm_ordinary: list of tuples (dim, (b,d)).
    :returns: list of tuples (dim, (b,d)) representing a persistence diagram.
    """
    dgm_relative = [(dim + 1, (-d if np.isfinite(d) else 0, -b)) for dim, (b, d) in dgm_ordinary]
    return dgm_relative


def plot_one_skeleton(X, x0, alpha):
    """Visualize the 1-skeleton of the VietorisRips complex built on `X`.
    """
    St = get_filtered_complex(X, x0, alpha, 1)

    simplices = [simplex for simplex, simplex_v in St.get_filtration() if (len(simplex) == 2)]
    n_points = X.shape[0]
    for simplex in simplices:
        i_a, i_b = simplex
        if (i_a < n_points) & (i_b < n_points):
            plt.plot([X[i_a][0], X[i_b][0]],
                     [X[i_a][1], X[i_b][1]], c='grey', alpha=0.3)
