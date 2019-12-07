#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This module provides facilities function in order to compute various quantities
related to the discrete curvature around a given point in R^3. The discrete
curvature is defined as the angular defect around a given vertex.


Let A be a given atom of a molecule and *(A) the set of atoms connected by
a bond to A. `coord` is a (N x 3) matrix each line being the cartesian
coordinates of one atom. The first line of the matrix corresponds to atom A
and the following lines corresponds to *(A) atoms.

Let us consider a set of point Reg*(A) corresponding to a set of
regularized point B for which the distances AB are equal. Let us consider
the plane P defined by the points Reg*(A). This plane is defined exactly
if Reg*(A) is an ensemble of 3 points.
The
pyramidalization angle is then define as `Pyr(A) = theta - pi / 2` where
theta is the angle


"""

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = ""
__email__ = "germain.vallverdu@univ-pau.fr"
__status__ = "Development"
__date__ = "June, 2018"

from .core import *
from .utils import *


