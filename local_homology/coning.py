import numpy as np

from sklearn.metrics import pairwise_distances

from gudhi import RipsComplex, SimplexTree


## Coning Variant
def is_point_in_ball(X, point, epsilon, return_distances=True):
    """Determine which points from X are in the epsilon-ball centered at point.
    :param X: coordinates of points.
    :type X: np.ndarray of size (n_pts, dim)
    :param point: center of the ball that defines the neighborhood, does not have to be in X.
    :type point: np.ndarray of shape (1, dim)
    :params epsilon: size of the neighborhood to consider
    :type epsilon: float, positive
    
    :returns: arrays specifying whether a given ball is, or not, in the complex.
    :
    """
    distance_from_point =  pairwise_distances(X, point, metric="euclidean")
    distance_to_boundary = epsilon - distance_from_point
    is_in_ball = (distance_to_boundary >= 0)[:,0]
    if return_distances:
        return is_in_ball, distance_to_boundary
    return is_in_ball

def build_local_complex(X, point, epsilon, dimensions):
    """Build a simplicial complex on vertices from
    X \cap B(point, \epsilon), along with an extra vertex representing the boundary.
    :param X: coordinates of points.
    :type X: np.ndarray of size (n_pts, dim)
    :param point: center of the ball that defines the neighborhood, does not have to be in X.
    :type point: np.ndarray of shape (1, dim)
    :params epsilon: size of the neighborhood to consider
    :type epsilon: float, positive
    :param dimensions: list of dimensions to compute the homology of."""
    X_dist = pairwise_distances(X, metric="euclidean")
    is_in_ball, distance_to_boundary = is_point_in_ball(X, point,
                                                        epsilon, return_distances=True)
    
    X_local = X[is_in_ball]
    X_dist_local = X_dist[is_in_ball][:, is_in_ball]
    local_boundary = distance_to_boundary[is_in_ball]

    X_local_stared = np.concatenate([X_dist_local, local_boundary.T], axis=0)
    X_local_stared = np.concatenate([X_local_stared,
                                     np.concatenate([local_boundary, np.array([[0]])], axis=0)], axis=1)
    
    rc = RipsComplex(distance_matrix=X_local_stared)
    st = rc.create_simplex_tree(max_dimension=max(dimensions)+1)
    n_points_in_ball = np.sum(is_in_ball)
    star_vertex = n_points_in_ball # indexing starts at 0! 
    st.assign_filtration([star_vertex], -1)
    return st

def get_persistence(simplexTree, dimensions):
    """Calculate the persistent homology of the abstract simplicial complex,
    filtering by positive values and dimensions.
    :param simplexTree: a simplcial complex, as returned by `build_local_complex`
    :type simplexTree: simplexTree
    :param dimensions: list of dimensions of persistent homology.
    :type dimensions: list of ints.
    
    :returns: persistence diagram in the Gudhi format.
    :rtype: list of tuples (dim, (birth, death)).
    """
    pd = simplexTree.persistence()
    return [p for p in pd if (p[0] in dimensions)&(p[1][0]>=0.)]

