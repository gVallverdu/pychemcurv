#!/usr/bin/env python
# coding: utf-8


"""
This module implements the `VertexAtom` and `Hybridization` classes in order
to compute the local curvature and the hybridation, respectively.
"""

import numpy as np

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"

__all__ = ["VertexAtom", "Hybridization"]


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
    r"""Given a set of N points in :math:`\mathbb{R}^3`, compute an orthonormal basis of 
    vectors, the first two belonging to the plane and the third one being normal 
    to the plane. In the particular case where N equal 3, there is an exact
    definition of the plane as the three points define an unique plan.

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


class VertexAtom:
    r"""
    This object represents an atom (or a point) associated to a vertex of the
    squeleton of a molecule. Hereafter are reminded some definitions. To a 
    complete overview, take a read at the following publications J Phys Chem ...

    We denote by A a given atom caracterized by its cartesian coordinates 
    corresponding to a vector in :math:`\mathbb{R}^3`. This atom A is bonded to
    one or several atoms B. The atoms B, bonded to atoms A belong to 
    :math:`\star(A)` and are caracterized by their cartesian coordinates defined
    as vectors in :math:`\mathbb{R}^3`. The geometrical
    object obtained by drawing a segment between bonded atoms is called the
    skeleton of the molecule and is the initial geometrical picture for a molecule.
    This class is defined from the cartesian coordinates of atom A and the atoms
    belonging to :math:`\star(A)`.

    More generally, the classes only considers points in :math:`\mathbb{R}^3`
    corresponding to the cartesian coordinates of the atoms.
    In consequence, the class can be used for all cases
    where a set of point in :math:`\mathbb{R}^3` is relevant.

    The following quantities are computed. For a complete definition, please
    read the reference cited above (JPC 2020). All the quantities need that 
    there are at least 3 atoms in :math:`\star(A)`.

    pyramidalization angle ``pyrA``
        The pyramidalization angle, **in degrees**. :math:`pyrA = \theta - \pi/2`
        where :math:`\theta` is the angle between the normal vector of the plane
        containing the atoms B of :math:`\star(A)` and a vector along a bond 
        between atom A and one B atom.

        An exact definition of pyrA needs that A is bonded to exactly 3 atoms in 
        order to be able to define a uniq plane that contains the atoms B
        belonging to :math:`\star(A)`. Nevertheless, pyrA is computed if
        more than 3 atoms are bonded to atom A by computing the best fitting plane
        of atoms belonging to :math:`\star(A)`.

    pyramidalization angle, ``pyrA_r``
        The pyramidalization angle **in radians**.

    improper angle, ``improper``
        The improper angle corresponding to the dihedral angle between the 
        planes defined by atoms (i, j, k) and (j, k, l), atom i being atom A and
        atoms j, k and l being atoms of :math:`\star(A)`. In consequence, the
        improper angle is defined only if there are 3 atoms in :math:`\star(A)`.

        The value of the improper angle is returned in radians.

    angular defect, ``angular_defect``
        The angluar defect is defined as :math:`2\pi - \sum_{F\in\star(A)} \alpha_F`
        where :math:`\alpha_F` are the angles at the vertex A of the faces 
        :math:`F\in\star(A)`. The angular defect is computed whatever the number
        of atoms in :math:`\star(A)`.

        The value of the angular defect is returned in radians.

    spherical curvature, ``spherical_curvature``
        The spherical curvature is computed as the radius of the osculating
        sphere of atoms A and atoms belonging to :math:`\star(A)`. The
        spherical curvature is computed as

        .. math::

            \kappa(A) = \frac{1}{\sqrt{\ell^2 + \dfrac{OA^2 - \ell^2}{4z_A^2}}}

        where O is the center of the circumbscribed circle of atoms in 
        :math:`\star(A)` ; A the vertex atom ; OA the distance between O and A ;
        :math:`\ell` the distance between O and atoms B of :math:`\star(A)` ; 
        :math:`z_A` the distance of atom A to the plane defined by 
        :math:`\star(A)`. The spherical curvature is defined only if there are 
        3 atoms in :math:`\star(A)`.

    pyramidalization distance ``pyr_distance``
        Distance of atom A to the plane define by :math:`\star(A)` or
        the best fitting plane of :math:`\star(A)`. 

        The value of the distance is in the same unit as the coordinates.

    If the number of atoms B in :math:`\star(A)` is not suitable to compute some
    properties, `np.nan` is returned.

    Note that the plane defined by atoms B belonging to :math:`\star(A)` is exactly 
    defined *only* in the case where there are three atoms B in :math:`\star(A)`. 
    In the case of pyrA, if there are more than 3 atoms in :math:`\star(A)`, the
    class use the best fitting plane considering all atoms in :math:`\star(A)` and 
    compute the geometrical quantities.
    """

    def __init__(self, a, star_a):
        r"""
        Args:
            a (np.ndarray): cartesian coordinates of atom A in :math:`\mathbb{R}^3`
            star_a (nd.array): (N x 3) cartesian coordinates of atoms B in :math:`\star(A)`
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
        r""" Coordinates of atoms B belonging to :math:`\star(A)` """
        return self._star_a

    @property
    def reg_star_a(self):
        r"""
        Regularized coordinates of points B in :math:`\star(A)` such as all 
        distances between A and points B are equal to unity. This corresponds to 
        :math:`Reg_{\epsilon}\star(A)` with :math:`\epsilon` = 1.
        """

        u = self._star_a - self._a
        norm = np.linalg.norm(u, axis=1)
        u /= norm[:, np.newaxis]

        return self._a + u

    @property
    def improper(self):
        r"""
        Compute the improper angle in randians between planes defined by atoms 
        (i, j, k) and (j, k, l). Atom A, is atom i and atoms j, k and l belong 
        to :math:`\star(A)`.

        ::

                     l
                     |
                     i
                    /  \
                  j     k

        This quantity is available only if the length of :math:`\star(A)` is 
        equal to 3.
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
        r"""
        Compute the distance of atom A to the plane define by :math:`\star(A)` or
        the best fitting plane of :math:`\star(A)`. The unit of the distance is the
        same as the unit of the coordinates of A and :math:`\star(A)`.
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
        """ Return the pyramidalization angle in radians. """

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
        """ Return the pyramidalization angle in degrees. """
        return np.degrees(self.pyrA_r)

    @property
    def angular_defect(self):
        r"""
        Compute the angular defect as a measure of the discrete curvature around 
        the vertex, point A.

        The calculation first looks for the best fitting plane of points 
        belonging to :math:`\star(A)` and sorts that points in order to compute 
        the angles between the edges connected to the vertex (A).
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
        r"""
        Compute the spherical curvature associated to the osculating sphere of
        points A and points B belonging to :math:`\star(A)`.
        Here, we assume that there is exactly 3 atoms B in :math:`\star(A)`.
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

    def as_dict(self, radians=True):
        """ 
        Return a dict version of all the properties that can be computed using
        this class. 

        Args:
            radians (bool): if True, angles are returned in radians (default)
        """
        data = {
            "pyrA": self.pyrA_r if radians else self.pyrA,
            "spherical_curvature": self.spherical_curvature,
            "angular_defect": self.angular_defect if radians else np.degrees(self.angular_defect),
            "improper": self.improper if radians else np.degrees(self.improper),
            "pyr_distance": self.pyr_distance
        }
        return data

    @staticmethod
    def from_pyramid(length, theta, n_star_A=3, radians=False, perturb=None):
        r"""Set up a VertexAtom from an ideal pyramidal structure.
        Build an ideal pyramidal geometry given the angle theta and randomize
        the positions by adding a noise of a given magnitude. The vertex of the 
        pyramid is the point A, and :math:`\star(A)`. are the points linked to 
        the vertex. The size of :math:`\star(A)`. is at least 3.

        :math:`\theta` is the angle between the normal vector of the plane defined
        from :math:`\star(A)` and the bonds between A and :math:`\star(A)`. 
        The pyramidalisation angle is defined from :math:`\theta` such as

        .. math::

            pyrA = \theta - \frac{\pi}{2}

        Args:
            length (float): the bond lenght
            theta (float): Angle to define the pyramid
            n_star_A (int): number of point bonded to A the vertex of the pyramid.
            radian (bool): True if theta is in radian (default False)
            perturb (float): Give the width of a normal distribution from which
                random numbers are choosen and added to the coordinates.

        Returns:
            A VertexAtom instance
        """
        r_theta = theta if radians else np.radians(theta)
        if n_star_A < 3:
            raise ValueError(
                "n_star_A = {} and must be greater than 3.".format(n_star_A))

        # build an ideal pyramid
        IB = length * np.sin(r_theta)
        step_angle = 2 * np.pi / n_star_A
        coords = [[0, 0, -length * np.cos(r_theta)]]
        coords += [[IB * np.cos(iat * step_angle),
                    IB * np.sin(iat * step_angle),
                    0] for iat in range(n_star_A)]
        coords = np.array(coords, dtype=np.float)

        # randomize positions
        if perturb:
            coords[1:, :] += np.random.normal(0, perturb, size=(n_star_A, 3))

        return VertexAtom(coords[0], coords[1:])

    def __str__(self):
        """ str representatio of the vertex atom """
        s = "pyrA: {:.4f} degrees\n".format(self.pyrA)
        s += "size of *(A): {}\n".format(len(self.star_a))
        s += "Atom A:\n{}\n".format(self.a)
        s += "Atoms B in *(A):\n{}\n".format(self.star_a)
        return s

    def __repr__(self):
        """ representation of the vertex atom """
        return "VertexAtom(a={}, star_a={})".format(self.a, self.star_a)

    def write_file(self, species="C", filename="vertex.xyz"):
        r"""Write the coordinates of atom A and atoms :math:`\star(A)`
        in a file in xyz format. You can set the name of species or a list but 
        the length of the list must be equal to the number of atoms.
        If filename is None, returns the string corresponding to the xyz file.

        Args:
            species (str, list): name of the species or list of the species names
            filename (str): path of the output file or None to get a string

        Returns:
            None if filename is a path, else, the string corresponding to the
            xyz file.
        """
        nat = len(self.star_a) + 1
        if len(species) != nat:
            species = nat * "C"

        lines = "%d\n" % nat
        lines += "xyz file from pychemcurv, "
        lines += "pyrA = %.4f degrees\n" % self.pyrA
        lines += "%2s %12.6f %12.6f %12.6f\n" % (species[0], 
                                                 self.a[0], self.a[1], self.a[2])
        for iat in range(0, nat - 1):
            lines += "%2s " % species[iat]
            lines += " ".join(["%12.6f" % x for x in self.star_a[iat]])
            lines += "\n"

        if filename is not None:
            with open(filename, "w", encoding="utf-8") as f:
                f.write(lines)
        else:
            return lines


class Hybridization:
    r"""
    This class compute the hybridization of the s and p atomic orbitals of a
    given atom A, considering the pyramidalization angle. All the properties
    are computed from the value of the pyramidalization angle.

    For a precise definition of the various quantities look at the definitions
    in JCP 2020.

    A chemical picture of the hybridization can be draw by considering the
    contribution of the :math:`p_z` atomic oribtal to the system :math:`\sigma`,
    or the contribution of the s atomic orbital to the system :math:`\pi`. 
    """

    def __init__(self, pyrA=None, vertex=None, radians=False):
        r"""
        Define the vertex atom or the pyramidalization angle of this atom.
        The pyramidalization angle `pyrA` value is considered first.
        If not provided, the pyramidalization angle is computed from the 
        definition of the vertex using the `VertexAtom` class.

        Args:
            pyrA (float): value of the pyramidalization angle
            vertex (VertexAtom): Vertex atom A and atoms of :math:`\star(A)`
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
        r""" 
        Value of :math:`c_{\pi}` 

        .. math::

            c_{\pi} = \sqrt{3} \tan Pyr(A)
        """
        return np.sqrt(2) * np.tan(self._pyrA)

    @property
    def lambda_pi(self):
        r""" 
        value of :math:`\lambda_{\pi}` 

        .. math::

            \lambda_{\pi} = \sqrt{1 - 3 \tan^2 Pyr (A)} 
        """

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
        r""" 
        value of hybridization number m 

        .. math::

            m = \left(\frac{c_{\pi}}{\lambda_{\pi}}\right)^2
        """
        return (self.c_pi / self.lambda_pi) ** 2

    @property
    def n(self):
        """ 
        value of hybridization number n 

        .. math::

            n = 3m + 2
        """
        return 3 * self.m + 2

    @property
    def hybridization(self):
        r""" 
        Compute the hybridization such as 

        .. math::

            s p^{(2 + c_{\pi}^2) / (1 - c_{\pi}^2)}

        This quantity corresponds to the amount of pz AO in the system 
        :math:`\sigma` and corresponds to the :math:`\tilde{n}` value defined 
        by Haddon.
        """
        return (2 + self.c_pi ** 2) / (1 - self.c_pi ** 2)

    def as_dict(self):
        r""" 
        Return a dict version of all the properties that can be computed with
        this class. Note that in the case of :math:`\lambda_{\pi}` and 
        :math:`c_{\pi}` the squared values are returned as as they are more 
        meaningfull.

        All quantities do not have unit.
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
