# coding: utf-8

"""
This module implements utility functions to compute several geometric 
properties.
"""

import numpy as np

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"


def center_of_mass(coords, masses=None):
    r"""Compute the center of mass of the points at coordinates `coords` with
    masses `masses`.

    Args:
        coords (np.ndarray): (N, 3) matrix of the points in :math:`\mathbb{R}^3`
        masses (np.ndarray): vector of length N with the masses

    Returns:
        The center of mass as a vector in :math:`\mathbb{R}^3`
    """
    # check coord array
    try:
        coords = np.array(coords, dtype=np.float64)
        coords = coords.reshape(coords.size // 3, 3)
    except ValueError:
        print("coords = ", coords)
        raise ValueError("Cannot convert coords in a numpy array of floats"
                         " with a shape (N, 3).")

    # check masses
    if masses is None:
        masses = np.ones(coords.shape[0])
    else:
        try:
            masses = np.array(masses, dtype=np.float64)
            masses = masses.reshape(coords.shape[0])
        except ValueError:
            print("masses = ", masses)
            raise ValueError("Cannot convert masses in a numpy array of "
                             "floats with length coords.shape[0].")
    if masses is None:
        masses = np.ones(coords.shape[0])

    return np.sum(coords * masses[:, np.newaxis], axis=0) / masses.sum()


def circum_center(coords):
    r"""Compute the coordinates of the center of the circumscribed circle from
    three points A, B and C in :math:`\mathbb{R}^3`.

    Args:
        coords (ndarray): (3x3) cartesian coordinates of points A, B and C. 

    Returns
        The coordinates of the center of the cicumscribed circle
    """
    try:
        coords = np.array(coords, dtype=np.float64).reshape(3, 3)
    except ValueError:
        print("coords = ", coords)
        raise ValueError("Cannot convert coords in a numpy array of floats"
                         " with a shape (3, 3).")

    # get coords of poins A, B and C
    a, b, c = coords

    # normal vector to ABC plane
    ABvAC = np.cross(b - a, c - a)

    # matrix M and vector B
    M = np.array([b - a, c - a, ABvAC])
    B = np.array([np.dot(b - a, (b + a) / 2),
                  np.dot(c - a, (c + a) / 2),
                  np.dot(ABvAC, a)])

    # solve linear system and return coordinates
    return np.dot(np.linalg.inv(M), B)


def get_plane(coords, masses=None):
    r"""Given a set of N points in :math:`\mathbb{R}^3`, compute an orthonormal 
    basis of vectors, the first two belonging to the plane and the third one 
    being normal to the plane. In the particular case where N equal 3, there is 
    an exact definition of the plane as the three points define an unique plan.

    If N = 3, use a gram-schmidt orthonormalization to compute the vectors. If
    N > 3, the orthonormal basis is obtained from SVD.

    Args:
        coords (np.ndarray): (N, 3) matrix of the points in :math:`\mathbb{R}^3`
        masses (np.ndarray): vector of length N with the masses

    Returns:
        Returns the orthonormal basis (vecx, vecy, n_a), vector n_a being 
        normal to the plane.
    """
    # check coord array
    try:
        coords = np.array(coords, dtype=np.float64)
        coords = coords.reshape(coords.size // 3, 3)
    except ValueError:
        print("coords = ", coords)
        raise ValueError("Cannot convert coords in a numpy array of floats"
                         " with a shape (N, 3).")

    # check masses
    if masses is None:
        masses = np.ones(coords.shape[0])
    else:
        try:
            masses = np.array(masses, dtype=np.float64)
            masses = masses.reshape(coords.shape[0])
        except ValueError:
            print("masses = ", masses)
            raise ValueError("Cannot convert masses in a numpy array of "
                             "floats with length coords.shape[0].")

    com = center_of_mass(coords, masses)

    if coords.shape == (3, 3):
        # the plane is exactly defined from 3 points
        vecx = coords[1] - coords[0]
        vecx /= np.linalg.norm(vecx)

        # vecy, orthonormal with vecx
        vecy = coords[2] - coords[0]
        vecy -= np.dot(vecy, vecx) * vecx
        vecy /= np.linalg.norm(vecy)

        # normal vector
        n_a = np.cross(vecx, vecy)

    else:
        # get the best fitting plane from SVD.
        _, _, (vecx, vecy, n_a) = np.linalg.svd(coords - com)

    return vecx, vecy, n_a
