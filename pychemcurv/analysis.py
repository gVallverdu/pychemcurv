# coding: utf-8

"""
This module implements the `CurvatureAnalyze` class to perform curvature 
analyses on molecular or periodic structures.
"""

import numpy as np
import pandas as pd

from pymatgen.core import Molecule, Structure
from pymatgen.core.bonds import obtain_all_bond_lengths
from .core import VertexAtom, TrivalentVertex, POAV1, POAV2


__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__email__ = "germain.vallverdu@univ-pau.fr"

__all__ = ["CurvatureAnalyzer"]


class CurvatureAnalyzer:
    """ This class provides helpful methods to analyze the local curvature
    on all atoms of a given structure. The structure is either a molecule or 
    a periodic structure. Once the structure is read, the class determines the
    connectivity of the structure in order to define all vertices. The 
    connectivity is defined on a distance criterion.
    """

    def __init__(self, structure, bond_tol=0.2, rcut=2.5, bond_order=None):
        """ The class needs a pymatgen.Structure or pymatgen.Molecule object as
        first argument. The other arguments are used to defined if two atoms are
        bonded or not.

        Args:
            structure (Structure, Molecule): A Structure or Molecule pymatgen 
                objects
            bond_tol (float): Tolerance used to determine if two atoms are 
                bonded. Look at `pymatgen.core.CovalentBond.is_bonded`.
            rcut (float): Cutoff distance in case the bond is not not known
            bond_order (dict): Not yet implemented

        """
        if isinstance(structure, (Molecule, Structure)):
            self.structure = structure
        else:
            raise TypeError("structure must a Molecule or Structure pymatgen"
                            " object. type(structure) is: " + str(type(structure)))

        self.bond_tol = bond_tol
        self.rcut = rcut
        self.bond_order = bond_order

        # compute distance matrix one time. You must call only one time
        # structure.distance_matrix to save computational time
        self._distance_matrix = self.structure.distance_matrix

        # look for bonds and set vertices
        self._vertices = []
        self._bonds = set()
        self._vertices_idx = []
        self._get_vertex()

        # fill a DataFrame with datas
        self._data = pd.DataFrame([])
        self._compute_data()

    @property
    def vertices(self):
        """ List of vertices associated to each atom of the molecule """
        return self._vertices

    @property
    def bonds(self):
        """ Set of tuples of bonded atom index """
        return self._bonds

    @property
    def vertices_idx(self):
        r""" List of tuples of the indexes of atoms in each vetex. The first
        index is atom A, the following are atoms of :math:`\star(A)`. """
        return self._vertices_idx

    @property
    def data(self):
        """ 
        Return a Data Frame that contains all the geometric and hybridization
        data. 
        """
        return self._data

    @property
    def distance_matrix(self):
        """ Returns the distance matrix between all atoms. For periodic 
        structures, this returns the nearest image distances. """
        return self._distance_matrix

    def _get_vertex(self):
        """ Look for all vertex defined as atoms bonded to at least
        3 neighbors and set up a list of VertexAtom object."""

        vertices = list()
        vertices_idx = list()
        bonds = set()
        errors = set()
        for isite, site_i in enumerate(self.structure):
            atom_A = site_i.coords

            star_a = list()
            vertex_idx = [isite]
            for jsite, site_j in enumerate(self.structure):
                if isite == jsite:
                    continue

                # check if i and j are bonded
                distance = self._distance_matrix[isite, jsite]
                bonded = False
                try:
                    # look for bond length database of pymatgen
                    # equivalent to CovalentBonds.is_bonded but avoid to compute
                    # two times the bond length
                    ref_distances = obtain_all_bond_lengths(site_i.specie,
                                                            site_j.specie)
                    # TODO: use ref_distances from a bond order
                    for rcut in ref_distances.values():
                        if distance < (1 + self.bond_tol) * rcut:
                            bonded = True

                except ValueError as e:
                    errors.add(str(e))
                    bonded = distance <= self.rcut

                # increment *(A) if i and j are bonded
                if bonded:
                    star_a.append(site_j.coords)
                    vertex_idx.append(jsite)
                    bonds.add(tuple(sorted([isite, jsite])))

            # set up VertexAtom objects
            if len(star_a) >= 3:
                vertices.append(VertexAtom(atom_A, star_a))
                vertices_idx.append(tuple(vertex_idx))
            else:
                vertices.append(None)
                vertices_idx.append(tuple(vertex_idx))

        self._vertices = vertices
        self._vertices_idx = vertices_idx
        self._bonds = bonds

        if errors:
            print("errors\n", "\n".join(errors))
            print("Default cutoff of {} was used for the above bond".format(self.rcut))

    def _compute_data(self):
        """ Compute geometric and hybridation data for all vertex in the
        structure and store them in a DataFrame. """

        data = list()
        nan_array = np.empty(3)
        nan_array.fill(np.nan)
        for vertex, vertex_idx in zip(self.vertices, self.vertices_idx):
            if vertex is None:
                vdata = {}
                vdata["n_star_A"] = 0
            else:
                if len(vertex.star_a) == 3:
                    vertex = TrivalentVertex(vertex.a, vertex.star_a)

                    try:
                        poav1 = POAV1(vertex=vertex)
                        poav2 = POAV2(vertex=vertex)
                    except ValueError as e:
                        print("Unable to compute all data.", 
                              vertex.as_dict(radians=False))
                        print(e)
                    vdata = {
                        **poav1.as_dict(
                            radians=False, include_vertex=True, list_obj=True),
                        **poav2.as_dict(radians=False, list_obj=True)
                    }
                else:
                    vdata = vertex.as_dict(radians=False, list_obj=True)

            ia = vertex_idx[0]
            vdata.update(atom_idx=ia, species=self.structure[ia].specie.symbol)

            distances = [self.distance_matrix[ia, j] for j in vertex_idx[1:]]
            vdata.update({"ave. neighb. dist.": np.mean(distances)})

            data.append(vdata)

        self._data = pd.DataFrame(data)

    @staticmethod
    def from_file(path, periodic=None):
        """ Returns a CurvatureAnalyze object from the structure at the 
        given path. This method relies on the file format supported with 
        pymatgen Molecule and Structure classes.

        Supported formats for periodic structure include CIF, POSCAR/CONTCAR, 
        CHGCAR, LOCPOT, vasprun.xml, CSSR, Netcdf and pymatgen’s JSON serialized 
        structures.

        Supported formats for molecule include include xyz, gaussian input 
        (gjf|g03|g09|com|inp), Gaussian output (.out|and pymatgen’s JSON 
        serialized molecules.

        Args:
            path (str): Path to the structure file
            periodic (bool): if True, assume that the file correspond to a
                periodic structure. Default is None. The method tries to read
                the file, first from the Molecule class and second from the
                Structure class of pymatgen.        
        """
        if periodic is None:
            # try to read as a molecule
            try:
                structure = Molecule.from_file(path)
            except ValueError as e1:
                print("Cannot read file as a molecule.")
                # Try to read as a periodic structure
                try:
                    structure = Structure.from_file(path)
                except ValueError as e2:
                    print("Cannot read file as a periodic structure.")
                    print("Try as a molecule, error:", e1)
                    print("Try as a structure, error:", e2)
                    raise ValueError(
                        "Unable to load structure from file '%s'" % path)

        elif periodic:
            # Structure object
            structure = Structure.from_file(path)

        else:
            # Molecule object
            structure = Molecule.from_file(path)

        print("Read structure, done.")

        return CurvatureAnalyzer(structure)

    def get_molecular_data(self):
        """ 
        Set up a model data dictionnary that contains species, coordinates and
        bonds of the structure. This dictionnary can be used as model data for
        further visulization in bio-dash.
        """

        # set up json file
        model_data = {"atoms": [], "bonds": []}

        # structure part
        for iat, site in enumerate(self.structure):
            name = "%s%d" % (site.specie.symbol, iat + 1)
            model_data["atoms"].append({"name": name,
                                        "serial": iat,
                                        "element": site.specie.symbol,
                                        "positions": site.coords.tolist()})

        # bonds part
        for bond in self.bonds:
            iat, jat = bond
            model_data["bonds"].append(
                {"atom1_index": iat, "atom2_index": jat}
            )

        return model_data