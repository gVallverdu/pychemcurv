# coding: utf-8

import setuptools

__version__ = ""
with open("pychemcurv/__init__.py", "r") as f_version:
    for line in f_version:
        if "__version__" in line:
            exec(line.strip())
            break

setuptools.setup(
    version=__version__,

    # find_packages()
    packages=setuptools.find_packages(exclude=["pychemcurv-data"]),
)
