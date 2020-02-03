#!/usr/bin/env python
# coding: utf-8


"""
This is the core module that provides the functions who computes the
geometrical descriptors associated to the curvature.

Notations
---------


A
*(A)
point B belongs to *(A)
I : barycenter of points in *(A)

"""

import numpy as np
from enum import Enum, unique

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__version__ = "2020.01.29"
__email__ = "germain.vallverdu@univ-pau.fr"
__status__ = "Development"

__all__ = ["VertexAtom", "Units", "Hybridization"]


@unique
class Units(Enum):
    # angstrom = "A"
    # meter = "m"
    # picometer = "pm"
    # nanometer = "nm"
    degrees = "degrees"
    radians = "radians"

    def __str__(self):
        return self.value


def center_of_mass(coords, masses=None):
    """Compute the center of mass of the points at coordinates `coords` with
    masses `masses`.

    Args:
        coords (np.ndarray): (N, 3) matrix of the points in R^3
        masses (np.ndarray): vector of length N with the masses

    Returns:
        The center of mass as a vector in R^3
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
    """Compute the coordinates of the center of the circumscribed circle from
    three points A, B and C in R^3.

    Args:
        coords (ndarray): (3x3) cartesian coordinates of points A, B and C. 
            One point per line.

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
    """Given a set of N points in 3D space, compute an orthonormal basis of 
    vectors, the first two belonging to the plane and the third one being normal 
    to the plane. In the particular case where N equal 3, there is an exact
    definition of the plane as the three points define an unique plan.

    If N = 3, use a gram-schmidt orthonormalization to compute the vectors. If
    N > 3, the orthonormal basis is obtained from SVD.

    Args:
        coords (np.ndarray): (N, 3) matrix of the points in R^3
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


class VertexAtom:
    """
    This is an object representing an atom associated to a vertex of the
    squeleton of a molecule. This atom A, is supposed to be bonded to atoms B
    belonging to *(A). The object provides properties that are caracteristic to
    the local geometry and the local curvature around atom A in the molecule.

    Note that the plane defined by atoms B belonging to *(A) is exactly 
    defined *only* in the case where there are three atoms B in *(A). 
    Nevertheless, if there are more than 3 atoms in *(A), the class use the
    best fitting plane considering all atoms in *(A) and compute the geometrical
    quantities.

    The pyramidaliation angle, the angular defect and the pyramidalization
    distance needs at least 3 atoms to be computed. On the contrary, the 
    improper angle and the spherical curvature are computed only if there
    are exactly three atoms B belonging to *(A). If the number of atoms B is
    wrong, `np.nan` is returned.

    Angles are return in radians. The pyramidalization angle is return in 
    degrees (`self.pyrA`) or in radians (`self.pyrA_r`).
    """

    def __init__(self, a, star_a):
        """
        Define atom A, one vertex of the squeleton of a given molecule and atoms
        B belonging to *(A). The unit of the angles computed in this class are 
        is radians or degrees depending on `ang_unit`. The distances computed 
        in this class are in the same unit as the input coordinates.

        Args:
            a (np.ndarray): cartesian coordinates of atom a in R^3 
            star_a (nd.array): (N x 3) cartesian coordinates of atoms in *(A)
        """
        try:
            self._a = np.array(a, dtype=np.float64).reshape(3)
        except ValueError:
            print("a = ", a)
            raise ValueError("Cannot convert a in a numpy array of floats.")

        try:
            self._star_a = np.array(star_a, dtype=np.float64)
            self._star_a = self._star_a.reshape(self._star_a.size // 3, 3)
        except ValueError:
            print("*A, star_a = ", star_a)
            raise ValueError("Cannot convert star_a in a numpy array of floats"
                             " with a shape (N, 3).")

    @property
    def a(self):
        """ Coordinates of atom A """
        return self._a

    @property
    def star_a(self):
        """ Coordinates of atoms belonging to *(A) """
        return self._star_a

    @property
    def reg_star_a(self):
        """
        Regularize the coordinates of the points in *(A) such as all distances 
        between points belonging to *(A) and A are equal to unity.
        """

        u = self._star_a - self._a
        norm = np.linalg.norm(u, axis=1)
        u /= norm[:, np.newaxis]

        return self._a + u

    @property
    def improper(self):
        """
        Compute the improper angle between planes defined by atoms (i, j, k) and
        (j, k, l). Atom A, is atom i and atoms i, j and k belong to *(A).

                     l
                     |
                     i
                    / \
                  j     k

        This quantity is available only if the length of *(A) is equal to 3.
        """

        # improper angle is defined only in the case of 3 atoms in *(A)
        if self._star_a.shape[0] != 3:
            return np.nan

        # get coords
        icoords = self._a
        jcoords, kcoords, lcoords = self._star_a

        # compute vectors
        vij = jcoords - icoords
        vjk = kcoords - jcoords
        vlk = kcoords - lcoords
        m = np.cross(vij, vjk)  # perpendicular to ijk
        n = np.cross(vlk, vjk)  # perpendicular to jkl

        # print("vij ", vij)
        # print("vjk ", vjk)
        # print("vlk ", vlk)
        # print("m   ", m)
        # print("n   ", n)

        # compute the angle
        theta = np.arctan2(np.dot(vij, n) * np.linalg.norm(vjk), np.dot(m, n))
        # theta2 = np.arccos(np.dot(m, n) / np.linalg.norm(m) / np.linalg.norm(n))
        # print(np.degrees(theta), np.degrees(theta2))

        return theta

    @property
    def pyr_distance(self):
        """
        Compute the distance of atom A to the plane define by *(A) or
        the best fitting plane of *(A). The unit of the distance is the same 
        as the unit of the coordinates of A and *(A).
        """
        # pyramidalization distance needs at least 3 atoms in *(A)
        if self._star_a.shape[0] < 3:
            return np.nan

        # compute normal vector of *(A)
        _, _, n_a = get_plane(self._star_a)
        com = center_of_mass(self._star_a)

        return np.abs(np.dot(self._a - com, n_a))

    @property
    def pyrA_r(self):
        """ Return the pyramidalization angle in radians.
        A pyramidal structure is assumed such as point A is the vertex of the 
        pyramid and other points belong to *(A).
        """

        # check input coords
        if self._star_a.shape[0] < 3:
            return np.nan

        # get the normal vector to the plane defined from *(A)
        _, _, n_a = get_plane(self.reg_star_a)

        # change the direction of n_a to be the same as IA.
        IA = self._a - center_of_mass(self.reg_star_a)
        n_a = -n_a if np.dot(IA, n_a) < 0 else n_a

        # compute pyrA
        v = self.reg_star_a[0] - self._a
        v /= np.linalg.norm(v)
        pyrA = np.arccos(np.dot(v, n_a)) - np.pi / 2

        return pyrA

    @property
    def pyrA(self):
        """ Return the pyramidalization angle in degrees.
        A pyramidal structure is assumed such as point A is the vertex of the 
        pyramid and other points belong to *(A).
        """
        return np.degrees(self.pyrA_r)

    @property
    def angular_defect(self):
        """
        Compute the angular defect as a measure of the discrete curvature around 
        the vertex, point A in R^3, and points B in R^3, bonded to point A, 
        belonging to *(A).

        The calculation first looks for the best fitting plane of points 
        belonging to *(A) and sorts that points in order to compute the angles 
        between the edges connected to the vertex (A).
        """

        # check and regularize coords
        if self._star_a.shape[0] < 3:
            return np.nan

        # get P the plane of *(A)
        vecx, vecy, _ = get_plane(self.reg_star_a)

        # compute all angles with vecx in order to sort atoms of *(A)
        com = center_of_mass(self.reg_star_a)
        u = self.reg_star_a - com
        norm = np.linalg.norm(u, axis=1)
        u /= norm[:, np.newaxis]
        cos = np.dot(u, vecx)
        angles = np.where(np.dot(u, vecy) > 0,
                          np.arccos(cos),
                          2 * np.pi - np.arccos(cos))

        # sort points according to angles
        idx = np.arange(angles.size)
        idx = idx[np.argsort(angles)]
        idx = np.append(idx, idx[0])

        # compute curvature
        ang_defect = 2 * np.pi
        for i, j in np.column_stack([idx[:-1], idx[1:]]):
            u = self.reg_star_a[i, :] - self._a
            u /= np.linalg.norm(u)

            v = self.reg_star_a[j, :] - self._a
            v /= np.linalg.norm(v)

            cos = np.dot(u, v)
            ang_defect -= np.arccos(cos)

        return ang_defect

    @property
    def spherical_curvature(self):
        """
        Compute the spherical curvature associated to the osculating sphere at
        a give point A of a molecule bonded to atoms B belonging to *(A).
        Here, we assume that there is 3 atoms B in *(A).
        """

        # check length of *(A)
        if self._star_a.shape[0] != 3:
            return np.nan

        # plane *(A)
        point_O = circum_center(self._star_a)
        _, _, n_a = get_plane(self._star_a)

        # needed length
        l = np.linalg.norm(self._star_a[0] - point_O)
        z_A = np.dot(self._a - point_O, n_a)
        OA = np.linalg.norm(self._a - point_O)

        # spherical curvature
        if np.isclose(z_A, 0, atol=0, rtol=1e-7):
            kappa = np.nan
        else:
            kappa = 1 / np.sqrt(l**2 + (OA**2 - l**2)**2 / (4 * z_A**2))

        return kappa

    def as_dict(self, unit=Units.radians):
        """ Return a dict version of all local quantities computed in this
        class. """
        radians = unit == Units.radians
        data = {
            "pyrA": self.pyrA_r if radians else self.pyrA,
            "spherical_curvature": self.spherical_curvature,
            "angular_defect": self.angular_defect if radians else np.degrees(self.angular_defect),
            "improper": self.improper if radians else np.degrees(self.improper),
            "pyr_distance": self.pyr_distance
        }
        return data


class Hybridization:
    """
    This class compute the hybridization of the s and p atomic orbitals of a
    given atom A, considering the pyramidalization angle.
    """

    def __init__(self, pyrA=None, vertex=None, radians=False):
        """
        Define the vertex atom A or the pyramidalization angle of this atom.
        In priority, the pyrA value is used. If not provided, the 
        pyramidalization angle is computed from the definition of the vertex.

        Args:
            pyrA (float): value of the pyramidalization angle
            vertex (VertexAtom): atom A, corresponding to a vertex of a molecule
            radians (bool): if true the value of pyrA is in radian (default False)
        """
        if pyrA is not None:
            self._pyrA = pyrA if radians else np.radians(pyrA)
        elif vertex is not None:
            if isinstance(vertex, VertexAtom):
                self._pyrA = vertex.pyrA_r
            else:
                raise TypeError("vertex must be of type VertexAtom")
        else:
            raise ValueError(
                "You have to provide either pyrA or a vertex atom.")

    @property
    def pyrA(self):
        """ Pyramidalization angle in degrees """
        return np.degrees(self._pyrA)

    @property
    def pyrA_r(self):
        """ Pyramidalization angle in radians """
        return self._pyrA

    @property
    def c_pi(self):
        """ Value of c_pi """
        return np.sqrt(2) * np.tan(self._pyrA)

    @property
    def lambda_pi(self):
        """ value of lambda_pi """

        # check domain definition of lambda_pi
        values = 1 - 2 * np.tan(self._pyrA) ** 2
        wrong = values < 0
        if np.all(wrong):
            raise ValueError("lambda_pi is not define. "
                             "pyrA (degrees) = {}".format(self.pyrA))
        elif np.any(wrong):
            values = np.where(values > 0, values, np.nan)

        return np.sqrt(values)

    @property
    def m(self):
        """ value of hybridization number m """
        return (self.c_pi / self.lambda_pi) ** 2

    @property
    def n(self):
        """ value of hybridization number n """
        return 3 * self.m + 2

    @property
    def hybridization(self):
        """ 
        Compute the hybridization such as s p^{(2 + C_pi^2) / (1 - C_pi^2)}.
        This quantity corresponds to the amount of pz AO in the system sigma
        and corresponds to the $\tilde{n}$ value defined by Haddon.
        """
        return (2 + self.c_pi ** 2) / (1 - self.c_pi ** 2)

    def as_dict(self):
        """ 
        Return a dict version of all local quantities computed in this class. 
        The square values of lambda_pi and c_pi coefficients are returned as 
        they are more meaningfull.
        """
        data = {
            "hybridization": self.hybridization,
            "n": self.n,
            "m": self.m,
            # "lambda_pi": self.lambda_pi,
            # "c_pi": self.c_pi,
            "c_pi^2": self.c_pi**2,
            "lambda_pi^2": self.lambda_pi**2,
        }
        return data
