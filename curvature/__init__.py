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


import numpy as np

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = ""
__email__ = "germain.vallverdu@univ-pau.fr"
__status__ = "Development"
__date__ = "June, 2018"

__all__ = ["get_improper", "get_pyr_angle", "get_molecule_curvature",
           "get_discrete_curvature", "get_pyr_distance"]


def get_improper(coords):
    """
    Assuming i is bonded to j, k and l, compute the improper angle between planes
    defined by (i, j, k) and (j, k, l).


                  l
                  |
                  i
                 / \
                j   k

    Args:
        coords (4 x 3 array): Coordinnates of point i, j, k and l. i must be the first

    Returns:
        improper angle (degrees)
    """

    coords = np.array(coords)
    if coords.shape != (4, 3):
        raise ValueError("The shape of the input coordinates must be (4, 3) "
                         "corresponding to 4 poinst in R^3.")

    icoords, jcoords, kcoords, lcoords = coords

    v1 = kcoords - lcoords
    v2 = jcoords - kcoords
    v3 = icoords - jcoords
    v23 = np.cross(v2, v3)  # perpendicular to j, k, l
    v12 = np.cross(v1, v2)  # perpendicular to i, j, k

    theta = np.arctan2(np.linalg.norm(v2) * np.dot(v1, v23), np.dot(v12, v23))

    return np.degrees(theta)


def get_pyr_distance(coords):
    """
    Compute the distance between atom A and the average plane defined from
    atoms *(A). `coords` is the (N x 3) matrix of the cartesian coordinates of
    atoms A and *(A). Cartesian coordinates of atom A is the first line
    `coords[0, :]` and cartesian coordinates of atoms *(A) are the next lines
    `coords[1:, :]`.

    Args:
        coords (N x 3 array): cartesian coordinnates in R^3.

    Returns:
        distance (float) in the same unit as the input coordinates.
    """

    coords = np.array(coords)
    if coords.shape[0] < 4:
        raise ValueError("Input coordinates must correspond to at least 4"
                         " points in R^3.")
    if coords.shape[1] != 3:
        raise ValueError("Input coordinates must correspond to point in R^3.")

    # barycenter of the points *(A)
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Compute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    return np.dot(coords[0, :] - G, u_norm)


def get_pyr_distance2(coords):
    """
    Compute the distance between atom A and the average plane defined from
    atoms *(A). `coords` is the (N x 3) matrix of the cartesian coordinates of
    atoms A and *(A). Cartesian coordinates of atom A is the first line
    `coords[0, :]` and cartesian coordinates of atoms *(A) are the next lines
    `coords[1:, :]`.

    Args:
        coords (N x 3 array): cartesian coordinnates in R^3.

    Returns:
        distance (float) in the same unit as the input coordinates.
    """

    coords = np.array(coords)
    if coords.shape[0] < 4:
        raise ValueError("Input coordinates must correspond to at least 4"
                         " points in R^3.")
    if coords.shape[1] != 3:
        raise ValueError("Input coordinates must correspond to point in R^3.")

    # barycenter of the points *(A)
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Compute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    return np.dot(coords[0, :] - G, u_norm)

def get_pyr_angle(coords):
    """
    Assuming the first point is linked to the followings, compute the pyramidalization
    following the definition of Haddon et al. The anlge of pyramidalization is
    computed as the angle between the vector normal to the (j, k, l) plane and the
    bonds between i and j, k and l.

    Args:
        coords (4 x 3 array): Coordinnates of the in R^3.

    Returns:
        list of float corresponding to the angles
    """

    coords = np.array(coords)
    if coords.shape != (4, 3):
        raise ValueError("The shape of the input coordinates must be (4, 3) "
                         "corresponding to 4 poinst in R^3.")

    # barycenter of the points
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Cpmpute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    # change the
    GI = coords[0, :] - G
    if np.dot(GI, u_norm) < 0:
        u_norm = - u_norm

    angles = list()
    for coord in coords[1:]:
        v = coord - coords[0, :]
        v /= np.linalg.norm(v)
        cos = np.dot(u_norm, v)
        angles.append(np.degrees(np.arccos(cos)))

    haddon = np.array(angles).mean() - 90.0

    return haddon, angles


def get_molecule_curvature(mol, rcut=2.0):
    """
    Compute the curvature on all atoms in the molecule. For each atom, the
    discrete curvature is computed as the angular defect. Each atom is considered
    as a vertex, and the edges are determined by atoms bonded to the vertex. The
    list of edges is build using rcut as a distance cut-off.

    Args:
        mol (mg.Molecule): pymatgen molecule object
        rcut (float): a distance cut-off to determine the edges around the vertex

    Returns
        a Molecule object with the following properties:
            'curvature': the curvatures
            'neighbors': the number of neighbors of each atom
            'distances': the average distances of neighbors
    """

    mol = mol.copy()
    curvatures = list()
    nneighbors = list()
    distances = list()
    for atom in mol:
        coords = atom.coords
        n = 0
        ave_distance = 0
        for neighbor, distance in mol.get_neighbors(atom, rcut):
            coords = np.row_stack((coords, neighbor.coords))
            n += 1
            ave_distance += distance

        distances.append(ave_distance / n)
        nneighbors.append(n)
        curvatures.append(get_discrete_curvature(coords))

    mol.add_site_property("curvature", curvatures)
    mol.add_site_property("neighbors", nneighbors)
    mol.add_site_property("distances", distances)

    return mol


def get_discrete_curvature(coords):
    """
    Compute the discrete curvature (angular defect) around a vertex of a set of
    points in R^3. The function works for any number of points greater than 4:
    one vertex and 3 points surrounding the vertex defining the edges. The needed
    data are the list of the coordinates of the points, the first one being the
    vertex's one.
    The function first looks for the best fit plane of the points surrounding the
    vertex and sorts that points in order to compute the angles between the edges
    connected to the vertex and the the angular defect.

    Args:
        coords (ndarray N x 3): Coordinates of the points. The first one is the vertex.

    Returns
        curvature (float)
    """
    coords = np.array(coords)

    if coords.shape[1] != 3:
        raise ValueError("3D coordinates are needed.")

    npts = coords.shape[0]
    if npts < 4:
        raise ValueError("Not enough points to compute curvature")

    # barycenter of the points
    G = coords[1:, :].sum(axis=0) / coords[1:, :].shape[0]

    # Cpmpute the best plane from SVD
    _, _, (vecx, vecy, u_norm) = np.linalg.svd(coords[1:, :] - G)

    # compute all angles with vectorx
    angles = list()
    for points in coords[1:, :]:
        u = points - G
        u /= np.linalg.norm(u)

        cos = np.dot(vecx, u)
        if np.dot(vecy, u) > 0:
            angles.append(np.degrees(np.arccos(cos)))
        else:
            angles.append(360 - np.degrees(np.arccos(cos)))

    # sort points according to angles
    idx = np.arange(1, npts)
    idx = idx[np.argsort(angles)]
    idx = np.append(idx, idx[0])

    # compute curvature
    curvature = 360
    for i, j in np.column_stack([idx[:-1], idx[1:]]):
        u = coords[i, :] - coords[0, :]
        u /= np.linalg.norm(u)

        v = coords[j, :] - coords[0, :]
        v /= np.linalg.norm(v)

        cos = np.dot(u, v)
        curvature -= np.degrees(np.arccos(cos))

    return curvature
