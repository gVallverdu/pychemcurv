# coding: utf-8

import setuptools

version = {}
with open("pychemcurv/version.py") as fp:
    exec(fp.read(), version)

setuptools.setup(
    version=version['__version__'],

    # find_packages()
    packages=setuptools.find_packages(exclude=["pychemcurv-data"]),
)
