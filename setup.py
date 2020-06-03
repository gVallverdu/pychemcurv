# coding: utf-8

import setuptools

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"
__version__ = "2020.6.3"

with open("README.rst", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pychemcurv",
    version=__version__,
    author=__author__,
    author_email=__email__,
    url="https://github.com/gVallverdu/pychemcurv",

    # A short description
    description="Discrete and local curvature applied to chemistry and chemical reactivity",

    # long description
    long_description=long_description,
    long_description_content_type="text/x-rst",

    # requirements
    install_requires=[
        "numpy", "pandas", "pymatgen", "matplotlib",
    ],
    # extra requirements
    extras_require={
        # for nglview visualization in jupyter notebook
        "viz": ["jupyter", "ase", "nglview"],
        # to run the dash app locally
        "app": ["dash", "dash-bio"],
    },

    # find_packages()
    packages=setuptools.find_packages(exclude=["pychemcurv-data"]),

    # 
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires='>=3.6',
)