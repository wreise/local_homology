import setuptools

setuptools.setup(
    name="local_homology",
    version="0.0.0.1",
    author="WReise",
    author_email="wojciech.reise@inria.fr",
    description="Approximation of (approximations of) local homology",
    long_description=" Bla ",
    url="https://github.com/wreise/local_homology",
    packages=setuptools.find_packages(),
    install_requires=["numpy>=1.19.4",
                      "scipy>=1.5.4",
                      "matplotlib>=3.3.3",
                      "gudhi>=3.3.0",
                      "scikit-learn>=0.23.2"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
