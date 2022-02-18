# Local homology

Currently, two ways to calculate "local (persistent) homology" of a point cloud `X` at a point `x_0` are implemented.
The first is akin to the alpha-filtration, while the second one is the r-filtration.

## Installation
It is advised to start from a clean environment and install directly from github:
```
conda create -n localHom python=3.8
conda activate localHom
python -m pip install git+https://github.com/wreise/local_homology.git
```
One can also clone (or download this repository).
After the first two steps described above, navigate to this directory and
`python -m pip install -e .`

## Try it out
See `Tutorial.ipynb`

# TODO:
- explain the expanding filtration drawbacks and advantages.
- establish guarantees for the coning process -> No guarantees, we are calculating something else. Derive guarantees for that?
- implement the true alpha-filtration.