# Local homology

This repository contains code for calculating local persistent homology on point clouds.
It is a python implementation relying on Gudhi, with ideas taken from [1]. A Tutorial is available.

## Installation
It is advised to start from a clean environment and install directly from github:
```
conda create -n localHom python=3.8
conda activate localHom
python -m pip install git+https://github.com/wreise/local_homology.git
```
One can also clone (or download this repository).
After the first two steps described above, navigate to this directory and
```
python -m pip install --pre --extra-index https://pypi.anaconda.org/scipy-wheels-nightly/simple -e .
```

## Try it out
See `Tutorial.ipynb`.
The nuances between different filtration variants are explained below.

## Mathematically
Currently, two ways to calculate "local (persistent) homology" of a point cloud `X` at a point `x_0` are implemented.

The first is akin to the alpha-filtration proposed in [1] and is implemented by `compute_local_homology_alpha`.
We fix a neighborhood size `epsilon` and we consider the filtration
$(X_alpha \cap B(x_0,\epsilon), X_\alpha\cap \partial B(x_0,\epsilon))$.
We approximate $X_\alpha$ with the Rips complex and relative homology is treated with coning, by defining at what filtration value each point
$x\in X\cap B(x_0,\epsilon)$ is connected to the boundary. There are three implemented strategies:
`distance_to_point_outside_ball` (default), `distance_to_boundary`, `distance_to_expanding_boundary`.
An illustration is provided in `implemented_variants.png`.

The second one is exactly the r-filtration from [1]. It is implemented following their suggestions, by calculating the normal homology of the sublevel sets
of a certain function and applying a duality theorem. 

Future work should include establishing guarantees (or counter-examples) for homology inference using the proposed alpha-filtration.

1. Skraba, P. & Wang, B. Approximating Local Homology from Samples. in Proceedings of the Twenty-Fifth Annual ACM-SIAM Symposium on Discrete Algorithms 174–192 (Society for Industrial and Applied Mathematics, 2014). doi:10.1137/1.9781611973402.13.
