# coding: utf-8

"""
Module ``pychemcur.core`` implements several classes in order to represents a vertex of 
a molecular squeleton and compute geometrical and chemical indicators related
to the local curvature around this vertex.

A complete and precise definition of all the quantities computed in the 
classes of this module can be found in article [JCP2020]_.

.. [JCP2020] Julia Sabalot-Cuzzubbo, Germain Salvato Vallverdu, Didier Bégué 
    and Jacky Cresson *Relating the molecular topology and local geometry: 
    Haddon’s pyramidalization angle and the Gaussian curvature*, J. Chem. Phys. 
    **152**, 244310 (2020). https://aip.scitation.org/doi/10.1063/5.0008368

.. [POAV2] Julia Sabalot-Cuzzubbo, N. Cresson, Germain Salvato Vallverdu, Didier
    Bégué and Jacky Cresson *Haddon’s POAV2 vs POAV theory for non-planar molecules*, 
    J. Chem. Phys. **159**, 174109 (2023). https://aip.scitation.org/doi/10.1063/5.0170800

"""

import numpy as np
from scipy.linalg import null_space
from .geometry import get_plane, circum_center, center_of_mass, get_dihedral

__all__ = ["VertexAtom", "TrivalentVertex", "POAV1", "POAV2"]


class VertexAtom:
    r"""
    This class represents an atom (or a point) associated to a vertex of the
    squeleton of a molecule. The used notations are the following. 
    We denote by A a given atom caracterized by its cartesian coordinates 
    corresponding to a vector in :math:`\mathbb{R}^3`. This atom A is bonded to
    one or several atoms B. The atoms B, bonded to atoms A belong to 
    :math:`\star(A)` and are caracterized by their cartesian coordinates defined
    as vectors in :math:`\mathbb{R}^3`. The geometrical
    object obtained by drawing a segment between bonded atoms is called the
    skeleton of the molecule and is the initial geometrical picture for a molecule.
    This class is defined from the cartesian coordinates of atom A and the atoms
    belonging to :math:`\star(A)`.

    More generally, the classes only considers points in :math:`\mathbb{R}^3`.
    The is not any chemical consideration here. In consequence, the class can be
    used for all cases where a set of point in :math:`\mathbb{R}^3` is relevant.
    """

    def __init__(self, a, star_a):
        r"""
        Args:
            a (np.ndarray): cartesian coordinates of point/atom A in :math:`\mathbb{R}^3`
            star_a (nd.array): (N x 3) cartesian coordinates of points/atoms B in :math:`\star(A)`
        """
        # check point/atom A
        try:
            self._a = np.array(a, dtype=np.float64).reshape(3)
        except ValueError:
            print("a = ", a)
            raise ValueError("Cannot convert a in a numpy array of floats.")

        # check points/atoms B in *(A)
        try:
            self._star_a = np.array(star_a, dtype=np.float64)
            self._star_a = self._star_a.reshape(self._star_a.size // 3, 3)
        except ValueError:
            print("*A, star_a = ", star_a)
            raise ValueError("Cannot convert star_a in a numpy array of floats"
                             " with a shape (N, 3).")

        if self._star_a.shape[0] < 3:
            print("*A, star_a = ", star_a)
            raise ValueError("The shape of *(A) is not relevant. Needs at least"
                             " 3 points/atoms in *(A)")

        # compute the regularized coordinates of atoms/points B in *(A)
        u = self._star_a - self._a
        self._distances = np.linalg.norm(u, axis=1)
        u /= self._distances[:, np.newaxis]
        self._reg_star_a = self._a + u

        # center of mass of atoms/points B in *(A)
        self._com = center_of_mass(self._star_a)

        # compute a normal vector of *(A)
        _, _, self._normal = get_plane(self._star_a)

        # compute a normal vector of the plane Reg *(A) using the regularized
        # coordinates of atoms/points B in *(A)
        _, _, self._reg_normal = get_plane(self._reg_star_a)

        # make the direction IA and the normal vectors of *(A) or Reg *(A) the same
        # I is the center of mass of *(A)
        IA = self.a - self.com
        if np.dot(IA, self._normal) < 0:
            self._normal = -self._normal
        if np.dot(IA, self.reg_normal) < 0:
            self._reg_normal = -self.reg_normal

    @staticmethod
    def from_pyramid(length, theta, n_star_A=3, radians=False, perturb=None):
        r"""Set up a VertexAtom from an ideal pyramidal structure.
        Build an ideal pyramidal geometry given the angle theta and randomize
        the positions by adding a noise of a given magnitude. The vertex of the 
        pyramid is the point A and :math:`\star(A)`. are the points linked to 
        the vertex. The size of :math:`\star(A)`. is at least 3.

        :math:`\theta` is the angle between the normal vector of the plane defined
        from :math:`\star(A)` and the bonds between A and :math:`\star(A)`. 
        The pyramidalisation angle is defined from :math:`\theta` such as

        .. math::

            pyrA = \theta - \frac{\pi}{2}

        Args:
            length (float): the bond length
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
        coords = np.array(coords, dtype=np.float64)

        # randomize positions
        if perturb:
            coords[1:, :] += np.random.normal(0, perturb, size=(n_star_A, 3))

        return VertexAtom(coords[0], coords[1:])

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
        Regularized coordinates of atoms/points B in :math:`\star(A)` such as all 
        distances between A and points B are equal to unity. This corresponds to 
        :math:`Reg_{\epsilon}\star(A)` with :math:`\epsilon` = 1.
        """
        return self._reg_star_a

    @property
    def normal(self):
        r"""
        Unitary vector normal to the plane or the best fitting plane of
        atoms/points Bi in :math:`\star(A)`.
        """
        return self._normal

    @property
    def reg_normal(self):
        r"""
        Unitary vector normal to the plane or the best fitting plane of
        atoms/points :math:`Reg B_i` in :math:`\star(A)`.
        """
        return self._reg_normal

    @property
    def com(self):
        r""" Center of mass of atoms/points B in :math:`\star(A)` """
        return self._com

    @property
    def distances(self):
        r"""
        Return all distances between atom A and atoms B belonging to 
        :math:`\star(A)`. Distances are in the same order as the atoms in 
        ``vertex.star_a``.
        """
        return self._distances

    def get_angles(self, radians=True):
        r"""
        Compute angles theta_ij between the bonds ABi and ABj, atoms Bi and
        Bj belonging to :math:`\star(A)`. The angle theta_ij is made by the 
        vectors ABi and ABj in the affine plane defined by this two vectors and 
        atom A. The computed angles are such as bond ABi are in a consecutive
        order.

        Args:
            radians (bool): if True (default) angles are returned in radians
        """
        if self._star_a.shape[0] == 3:
            angles = dict()
            for i, j in [(0, 1), (0, 2), (1, 2)]:
                u = self.reg_star_a[i, :] - self._a
                v = self.reg_star_a[j, :] - self._a

                cos = np.dot(u, v)
                if radians:
                    angles[(i, j)] = np.arccos(cos)
                else:
                    angles[(i, j)] = np.degrees(np.arccos(cos))

        else:
            # get P the plane of *(A)
            vecx, vecy, _ = get_plane(self.reg_star_a)

            # compute all angles with vecx in order to sort atoms of *(A)
            com = center_of_mass(self.reg_star_a)
            u = self.reg_star_a - com
            norm = np.linalg.norm(u, axis=1)
            u /= norm[:, np.newaxis]
            cos = np.dot(u, vecx)
            angles = np.where(np.dot(u, vecy) > 0, np.arccos(cos),
                              2 * np.pi - np.arccos(cos))

            # sort points according to angles
            idx = np.arange(angles.size)
            idx = idx[np.argsort(angles)]
            idx = np.append(idx, idx[0])

            # compute curvature
            angles = dict()
            for i, j in np.column_stack([idx[:-1], idx[1:]]):
                u = self.reg_star_a[i, :] - self._a
                u /= np.linalg.norm(u)

                v = self.reg_star_a[j, :] - self._a
                v /= np.linalg.norm(v)

                cos = np.dot(u, v)
                if radians:
                    angles[(i, j)] = np.arccos(cos)
                else:
                    angles[(i, j)] = np.degrees(np.arccos(cos))

        return angles

    @property
    def angular_defect(self):
        r"""
        Compute the angular defect in radians as a measure of the discrete 
        curvature around the vertex, point A.

        The calculation first looks for the best fitting plane of points 
        belonging to :math:`\star(A)` and sorts that points in order to compute 
        the angles between the edges connected to the vertex (A). See the
        get_angles method.
        """
        angles = self.get_angles(radians=True)
        ang_defect = 2 * np.pi - sum(angles.values())

        return ang_defect

    @property
    def pyr_distance(self):
        r"""
        Compute the distance of atom A to the plane define by :math:`\star(A)` or
        the best fitting plane of :math:`\star(A)`. The unit of the distance is the
        same as the unit of the coordinates of A and :math:`\star(A)`.
        """
        return np.abs(np.dot(self._a - self.com, self.normal))

    def as_dict(self, radians=True, list_obj=False):
        """ 
        Return a dict version of all the properties that can be computed using
        this class. Use `list_obj=True` to get a valid JSON object.

        Args:
            radians (bool): if True, angles are returned in radians (default)
            list_obj (bool): if True, numpy arrays are converted into list object (default False)

        Returns:
            A dict
        """
        angles = self.get_angles(radians=radians)
        angles = [[i, j, angle] for (i,j), angle in angles.items()]
        data = {
            "atom_A": self.a.tolist() if list_obj else self.a,
            "star_A": self.star_a.tolist() if list_obj else self.star_a,
            "reg_star_A": self.reg_star_a.tolist() if list_obj else self.reg_star_a,
            "distances": self.distances.tolist() if list_obj else self.distances,
            "angles": angles,
            "n_star_A": len(self.star_a),
            "angular_defect": self.angular_defect if radians else np.degrees(self.angular_defect),
            "pyr_distance": self.pyr_distance,
        }
        return data

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
        lines += "xyz file from pychemcurv\n"
        lines += "%2s %12.6f %12.6f %12.6f\n" % (species[0],
                                                 self.a[0], self.a[1], self.a[2])
        for iat in range(1, nat):
            lines += "%2s " % species[iat]
            lines += " ".join(["%12.6f" % x for x in self.star_a[iat - 1]])
            lines += "\n"

        if filename is not None:
            with open(filename, "w", encoding="utf-8") as f:
                f.write(lines)
        else:
            return lines

    def __str__(self):
        """ str representatio of the vertex atom """
        s = "angular defect: {:.4f} degrees\n".format(
            np.degrees(self.angular_defect))
        s += "size of *(A): {}\n".format(len(self.star_a))
        s += "Atom A:\n{}\n".format(self.a)
        s += "Atoms B in *(A):\n{}\n".format(self.star_a)
        return s

    def __repr__(self):
        """ representation of the vertex atom """
        return "VertexAtom(a={}, star_a={})".format(self.a, self.star_a)


class TrivalentVertex(VertexAtom):
    r"""
    This object represents an atom (or a point) associated to a vertex of the
    squeleton of a molecule bonded to exactly 3 other atoms (or linked to 3 
    other points). This correspond to the trivalent case.

    We denote by A a given atom caracterized by its cartesian coordinates 
    corresponding to a vector in :math:`\mathbb{R}^3`. This atom A is bonded to
    3 atoms B. The atoms B, bonded to atom A belong to 
    :math:`\star(A)` and are caracterized by their cartesian coordinates defined
    as vectors in :math:`\mathbb{R}^3`. The geometrical
    object obtained by drawing a segment between bonded atoms is called the
    skeleton of the molecule and is the initial geometrical picture for a molecule.
    This class is defined from the cartesian coordinates of atom A and the atoms
    belonging to :math:`\star(A)`.

    More generally, the classes only considers points in :math:`\mathbb{R}^3`.
    The is not any chemical consideration here. In consequence, the class can be
    used for all cases where a set of point in :math:`\mathbb{R}^3` is relevant.

    The following quantities are computed according the reference [JCP2020]_

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
        The angluar defect is defined as 

        .. math:

            2\pi - \sum_{F\in\star(A)} \alpha_F

        where :math:`\alpha_F` are the angles at the vertex A of the faces 
        :math:`F\in\star(A)`. The angular defect is computed whatever the number
        of atoms in :math:`\star(A)`.

        The value of the angular defect is returned in radians.

    spherical curvature, ``spherical_curvature``
        The spherical curvature is computed as the radius of the osculating
        sphere of atoms A and atoms belonging to :math:`\star(A)`. The
        spherical curvature is computed as

        .. math::

            \kappa(A) = \frac{1}{\sqrt{\ell^2 + \dfrac{(OA^2 - \ell^2)^2}{4z_A^2}}}

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
            a (np.ndarray): cartesian coordinates of point/atom A in :math:`\mathbb{R}^3`
            star_a (nd.array): (N x 3) cartesian coordinates of points/atoms B in :math:`\star(A)`
        """
        super().__init__(a, star_a)

        if self._star_a.shape[0] != 3:
            raise ValueError("The number of atoms/points in *(A) must be 3."
                             " star_a.shape is {}".format(self._star_a.shape))

    @staticmethod
    def from_pyramid(length, theta, radians=False, perturb=None):
        r"""Set up a VertexAtom from an ideal pyramidal structure.
        Build an ideal pyramidal geometry given the angle theta and randomize
        the positions by adding a noise of a given magnitude. The vertex of the 
        pyramid is the point A and :math:`\star(A)`. are the points linked to 
        the vertex. The size of :math:`\star(A)`. is 3.

        :math:`\theta` is the angle between the normal vector of the plane defined
        from :math:`\star(A)` and the bonds between A and :math:`\star(A)`. 
        The pyramidalisation angle is defined from :math:`\theta` such as

        .. math::

            pyrA = \theta - \frac{\pi}{2}

        Args:
            length (float): the bond length
            theta (float): Angle to define the pyramid
            radian (bool): True if theta is in radian (default False)
            perturb (float): Give the width of a normal distribution from which
                random numbers are choosen and added to the coordinates.

        Returns:
            A TrivalentVertex instance
        """
        va = VertexAtom.from_pyramid(
            length, theta, n_star_A=3, radians=radians, perturb=perturb
        )
        return TrivalentVertex(a=va.a, star_a=va.star_a)

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
        return get_dihedral(np.concatenate((self._a[np.newaxis, :], self._star_a)))

    @property
    def pyrA_r(self):
        """ Return the pyramidalization angle in radians. """

        # compute pyrA
        v = self.reg_star_a[0] - self._a
        v /= np.linalg.norm(v)
        pyrA = np.arccos(np.dot(v, self.reg_normal)) - np.pi / 2

        return pyrA

    @property
    def pyrA(self):
        """ Return the pyramidalization angle in degrees. """
        return np.degrees(self.pyrA_r)

    @property
    def spherical_curvature(self):
        r"""
        Compute the spherical curvature associated to the osculating sphere of
        points A and points B belonging to :math:`\star(A)`.
        Here, we assume that there is exactly 3 atoms B in :math:`\star(A)`.
        """

        # plane *(A)
        point_O = circum_center(self._star_a)

        # needed length
        l = np.linalg.norm(self._star_a[0] - point_O)
        z_A = np.dot(self._a - point_O, self.normal)
        OA = np.linalg.norm(self._a - point_O)

        # spherical curvature
        if np.isclose(z_A, 0, atol=0, rtol=1e-7):
            kappa = np.nan
        else:
            kappa = 1 / np.sqrt(l**2 + (OA**2 - l**2)**2 / (4 * z_A**2))

        return kappa

    def as_dict(self, radians=True, list_obj=False):
        """ 
        Return a dict version of all the properties that can be computed using
        this class. Use `list_obj=True` to get a valid JSON object.  

        Args:
            radians (bool): if True, angles are returned in radians (default)
            list_obj (bool): if True, numpy arrays are converted into list object (default False)

        Returns:
            A dict.
        """
        data = super().as_dict(radians=radians, list_obj=list_obj)
        data.update({
            "pyrA": self.pyrA_r if radians else self.pyrA,
            "spherical_curvature": self.spherical_curvature,
            "improper": self.improper if radians else np.degrees(self.improper),
        })
        return data

    def __str__(self):
        """ str representation of the vertex atom """
        s = "pyrA: {:.4f} degrees\n".format(self.pyrA)
        s += "Atom A:\n{}\n".format(self.a)
        s += "Atoms B in *(A):\n{}\n".format(self.star_a)
        return s

    def __repr__(self):
        """ representation of the vertex atom """
        return "TrivalentVertex(a={}, star_a={})".format(self.a, self.star_a)


class POAV1:
    r"""
    In the case of the POAV1 theory
    the POAV vector has the property to make a constant angle with each bond
    connected to atom A.

    This class computes indicators related to the POAV1 theory of R.C. Haddon
    following the link established between pyrA and the hybridization of a
    trivalent atom in reference [JCP2020]_.

    A chemical picture of the hybridization can be drawn by considering the
    contribution of the :math:`p` atomic oribtals to the system :math:`\sigma`,
    or the contribution of the s atomic orbital to the system :math:`\pi`. This
    is achieved using the m and n quantities. For consistency with POAV2 class,
    the attributes, ``hybridization``, ``sigma_hyb_nbr`` and ``pi_hyb_nbr`` 
    are also implemented but return the same values.
    """

    def __init__(self, vertex):
        r"""
        POAV1 is defined from the local geometry of an atom at a vertex of the
        molecule's squeleton.

        Args:
            vertex (TrivalentVertex): the trivalent vertex atom
        """
        if isinstance(vertex, TrivalentVertex):
            self.vertex = vertex
        elif isinstance(vertex, VertexAtom):
            self.vertex = TrivalentVertex(vertex.a, vertex.star_a)
        else:
            raise TypeError("vertex must be of type VertexAtom or of type"
                            " TrivalentVertex. vertex is {}".format(type(vertex)))

    @property
    def pyrA(self):
        """ Pyramidalization angle in degrees """
        return self.vertex.pyrA

    @property
    def pyrA_r(self):
        """ Pyramidalization angle in radians """
        return self.vertex.pyrA_r

    @property
    def poav(self):
        """ Return a unitary vector along the POAV vector """
        return self.vertex.reg_normal

    @property
    def c_pi(self):
        r""" 
        Value of :math:`c_{\pi}` in the ideal case of a :math:`C_{3v}` 
        geometry. Equation (22), with :math:`c_{1,2} = \sqrt{2/3}`.

        .. math::

            c_{\pi} = \sqrt{2} \tan Pyr(A)
        """
        return np.sqrt(2) * np.tan(self.pyrA_r)

    @property
    def lambda_pi(self):
        r""" 
        value of :math:`\lambda_{\pi}` in the ideal case of a :math:`C_{3v}` 
        geometry. Equation (23), with :math:`c^2_{1,2} = 2/3`.

        .. math::

            \lambda_{\pi} = \sqrt{1 - 2 \tan^2 Pyr (A)} 
        """

        # check domain definition of lambda_pi
        value = 1 - 2 * np.tan(self.pyrA_r) ** 2
        if value < 0:
            raise ValueError("lambda_pi is not define. "
                             "pyrA (degrees) = {}".format(self.pyrA))
        else:
            return np.sqrt(value)

    @property
    def m(self):
        r""" 
        value of hybridization number m, see equation (44) 

        .. math::

            m = \left(\frac{c_{\pi}}{\lambda_{\pi}}\right)^2
        """
        return (self.c_pi / self.lambda_pi) ** 2

    @property
    def n(self):
        """ 
        value of hybridization number n, see equation (47)

        .. math::

            n = 3m + 2
        """
        return 3 * self.m + 2

    @property
    def pi_hyb_nbr(self):
        r""" This quantity measure the weight of the s atomic orbital with
        respect to the p atomic orbital in the :math:`h_{\pi}` hybrid orbital 
        along the POAV vector.

        This is equal to m.
        """
        return self.m

    @property
    def sigma_hyb_nbr(self):
        """ This quantity measure the weight of the p atomic orbitals with
        respect to s in the hi hybrid orbitals along the bonds with atom A.

        This is equal to n
        """
        return self.n

    @property
    def hybridization(self):
        r""" Compute the hybridization such as 

        .. math::

            s p^{(2 + c_{\pi}^2) / (1 - c_{\pi}^2)}

        This quantity corresponds to the amount of p AO in the system 
        :math:`\sigma`. This is equal to n and corresponds to the 
        :math:`\tilde{n}` value defined by Haddon.

        TODO: verifier si cette quantité est égale à n uniquement dans 
        le cas C3v.
        """
#        return self.n
        return (2 + self.c_pi ** 2) / (1 - self.c_pi ** 2)

    def as_dict(self, radians=True, include_vertex=False, list_obj=False):
        r""" Return a dict version of all the properties that can be 
        computed with this class. Note that in the case of 
        :math:`\lambda_{\pi}` and :math:`c_{\pi}` the squared values are
        returned as they are more meaningfull. Use `list_obj= True` to
        obtain a valid JSON object.

        Args:
            radians (bool): if True, angles are returned in radians (default)
            include_vertex (bool): if True, include also vertex data
            list_obj (bool): if True, numpy arrays are converted into list object (default False)

        Returns:
            A dict.
        """
        data = {
            "hybridization": self.hybridization,
            "n": self.n,
            "m": self.m,
            # "lambda_pi": self.lambda_pi,
            # "c_pi": self.c_pi,
            "c_pi^2": self.c_pi ** 2,
            "lambda_pi^2": self.lambda_pi ** 2,
            "poav": self.poav.tolist() if list_obj else self.poav,
        }
        if include_vertex:
            data.update(self.vertex.as_dict(
                radians=radians, list_obj=list_obj))

        return data


class POAV2:
    r""" In the case of the POAV2 theory the POAV2 vector on atom A is 
    such as the set of hybrid molecular orbitals :math:`{h_{\pi}, h_1, h_2, h_3}` 
    is orthogonal ; where the orbitals :math:`h_i` are hybrid orbitals 
    along the bonds with atoms linked to atom A and :math:`h_{\pi}` is 
    the orbital along the POAV2 :math:`\vec{u}_{\pi}` vector.

    This class computes indicators related to the POAV2 theory of 
    R.C. Haddon following the demonstrations in the reference [POAV2]_.
    """

    def __init__(self, vertex):
        r""" POAV1 is defined from the local geometry of an atom at a 
        vertex of the molecule's squeleton.

        Args:
            vertex (TrivalentVertex): the trivalent vertex atom
        """
        if isinstance(vertex, TrivalentVertex):
            self.vertex = vertex
        elif isinstance(vertex, VertexAtom):
            self.vertex = TrivalentVertex(vertex.a, vertex.star_a)
        else:
            raise TypeError("vertex must be of type VertexAtom or of type"
                            " TrivalentVertex. vertex is {}".format(type(vertex)))

        self.angles = self.vertex.get_angles(radians=True)

    @property
    def matrix(self):
        """ Compute and return the sigma-orbital hybridization numbers 
        n1, n2 and n3 """
        cos_01 = np.cos(self.angles[(0, 1)])
        cos_02 = np.cos(self.angles[(0, 2)])
        cos_12 = np.cos(self.angles[(1, 2)])

        ui = self.vertex.reg_star_a - self.vertex.a

        M = np.array([
            [ui[2, 0] * cos_01 - ui[1, 0] * cos_02,
             ui[2, 1] * cos_01 - ui[1, 1] * cos_02,
             ui[2, 2] * cos_01 - ui[1, 2] * cos_02],
            [ui[0, 0] * cos_12 - ui[2, 0] * cos_01,
             ui[0, 1] * cos_12 - ui[2, 1] * cos_01,
             ui[0, 2] * cos_12 - ui[2, 2] * cos_01],
            [ui[1, 0] * cos_02 - ui[0, 0] * cos_12,
             ui[1, 1] * cos_02 - ui[0, 1] * cos_12,
             ui[1, 2] * cos_02 - ui[0, 2] * cos_12]
        ])

        return M

    @property
    def u_pi(self):
        r"""
        Return vector :math:`u_{\pi}` as the basis of the zero space of the 
        matrix M. This unitary vector support the POAV2 vector.
        """
        u = null_space(self.matrix)
        rank = u.shape[1]
        if rank != 1:
            raise ValueError("The rank of the null space is not equal to 1. "
                             "The POAV2 u_pi vector may not exist. "
                             "rank = %d" % rank)

        u = u.ravel()

        # make the direction of u_pi the same as IA (and thus reg_normal)
        # I is the center of mass of *(A)
        IA = self.vertex.a - self.vertex.com
        if np.dot(IA, u) < 0:
            u *= -1
        return u

    @property
    def sigma_hyb_nbrs(self):
        r"""
        Compute and return the sigma-orbital hybridization numbers n1, n2 and n3.
        These quantities measure the weight of the p atomic orbitals with
        respect to s in each of the :math:`h_i` hybrid orbitals along the bonds 
        with atom A.
        """
        cos_01 = np.cos(self.angles[(0, 1)])
        cos_02 = np.cos(self.angles[(0, 2)])
        cos_12 = np.cos(self.angles[(1, 2)])

        n1 = - cos_12 / cos_01 / cos_02
        n2 = - cos_02 / cos_12 / cos_01
        n3 = - cos_01 / cos_02 / cos_12

        return n1, n2, n3

    @property
    def pi_hyb_nbr(self):
        r"""
        This quantity measure the weight of the s atomic orbital with
        respect to the p atomic orbital in the :math:`h_{\pi}` hybrid orbital
        along the POAV2 vector.
        """
        n = self.sigma_hyb_nbrs

        w_sigma = sum([1 / (1 + ni) for ni in n])
        m = 1 / w_sigma - 1

        return m

    @property
    def pyrA_r(self):
        r"""
        Compute the angles between vector :math:`u_{\pi}` and all the bonds 
        between atom A and atoms B in :math:`\star(A)`. 
        """
        ui = self.vertex.reg_star_a - self.vertex.a
        scal = np.dot(ui, self.u_pi)
        return np.arccos(scal)

    @property
    def pyrA(self):
        return np.degrees(self.pyrA_r)

    def as_dict(self, radians=True, include_vertex=False, list_obj=False):
        r""" 
        Return a dict version of all the properties that can be computed with
        this class. Use `list_obj= True` to obtain a valid JSON object.

        Args:
            radians (bool): if True, angles are returned in radians (default)
            include_vertex (bool): if True, include also vertex data
            list_obj (bool): if True, numpy arrays are converted into list object (default False)

        Returns:
            A dict.
        """
        data = {
            "pi_hyb_nbr": self.pi_hyb_nbr,
            "u_pi": self.u_pi.tolist() if list_obj else self.u_pi,
            "matrix": self.matrix.tolist() if list_obj else self.matrix,
        }
        data.update({"n_%d" % i: ni
                     for i, ni in enumerate(self.sigma_hyb_nbrs, 1)})

        if include_vertex:
            data.update(self.vertex.as_dict(
                radians=radians, list_obj=list_obj))

        return data
