# coding: utf-8

import warnings
import numpy as np
import pychemcurv

import pandas as pd

__author__ = "Germain Salvato Vallverdu"
__copyright__ = ""
__version__ = "0.1"
__maintainer__ = "Germain Salvato Vallverdu"
__email__ = "germain.vallverdu@univ-pau.fr"

BOND_CUTOFFS = {
    # distances in angstrom
    ("H", "H"): 1.0,
    ("C", "C"): 1.7,
    ("C", "O"): 1.7,
    ("C", "S"): 1.9,
    ("C", "N"): 1.7,
    ("C", "Cl"): 1.9,
    ("C", "Br"): 2.2,
    ("C", "F"): 1.7,
    ("C", "H"): 1.2,
    ("C", "P"): 1.9,
    ("C", "H"): 1.2,
    ("N", "H"): 1.2,
    ("O", "H"): 1.2,
    ("S", "H"): 1.3,
    ("P", "H"): 1.3,
    ("P", "O"): 1.8,
    ("P", "N"): 1.8,
    ("S", "O"): 1.7,
    ("F", "H"): 1.2,
}

def get_bond_cutoff(specie1, specie2, rcut=2.5):
    """
    Return the bond cut-off in angstrom corresponding to the BOND_CUTOFF_DICT
    """
    if (specie1, specie2) in BOND_CUTOFFS:
        return BOND_CUTOFFS[(specie1, specie2)]
    elif (specie2, specie1) in BOND_CUTOFFS:
        return BOND_CUTOFFS[(specie2, specie1)]
    else:
        warnings.warn(f"bond {specie1}-{specie2} not known", category=UserWarning)
        return rcut


def read_xyz(file):
    """
    Read an xyz like file.
    The file is suposed to display the number of atom on the first line,
    followed by a title line and followed by the structure in cartesian
    coordinates. Each line contains the element as first item and the
    cartesian coordinates as 2d, 3th and 4th items. Example:

        3
        H2O molecule
        O   -0.111056  0.033897  0.043165
        H    0.966057  0.959148 -1.089095
        H    0.796629 -1.497157  0.403985

    If additional data are provided on each line, they are returned as
    atomic properties. Example, this file:

        3
        H2O molecule
        O   -0.111056  0.033897  0.043165   -1.8
        H    0.966057  0.959148 -1.089095    0.9
        H    0.796629 -1.497157  0.403985    0.9

    will return an `atomic_prop` dict such as:

        {"prop1": [-1.8, 0.9, 0.9]}

    Args:
        file (IO, string): An IO object with a readline() method or string 
            corresponding to a file path.

    Returns:
        * species, the list of elements
        * coords, a numpy array of floats with the shape (natom, 3)
        * None or atomic_prop dict such as {"prop1": [...], "prop2": [...]}
    """

    if isinstance(file, str):
        file = open(file, "r")

    # natom
    natom = int(file.readline().split()[0])
    _ = file.readline()

    # look for addition atomic properties
    nval_per_lines = [len(file.readline().split()) for _ in range(natom)]
    nprop = max(nval_per_lines)
    if all([nprop == nval for nval in nval_per_lines]):
        atomic_prop = list()
        nprop = nprop - 4
    else:
        atomic_prop = None
        nprop = 0

    # read the file
    file.seek(0)
    file.readline()  # natom
    file.readline()  # title

    coords = list()
    species = list()
    for _ in range(natom):
        data = file.readline().split()
        species.append(data[0])
        coords.append(data[1:4])

        if nprop > 0:
            atomic_prop.append(data[4:4 + nprop])

    # close file
    file.close()

    if nprop > 0:
        atomic_prop = np.array(atomic_prop, dtype=np.float).transpose()
        atomic_prop = {"prop%d" % i: d for i, d in enumerate(atomic_prop)}

    return species, np.array(coords, dtype=np.float), atomic_prop


def compute_data(species, coords, rcut=2.5, distances=None):
    """
    Compute the data over the whole molecule. For each atom, the following data 
    are computed:
    * pyramidalization angle, Pyr(A) (degrees)
    * angular defect (degrees)
    * spherical curvature, kappa, (A^-1)
    * improper angle (degrees)
    * distance from the average plane of neighbors
    * the average distances with neighbors
    * the number of neighbors
    * hybridization coefficients c_pi^2 and lambda_pi^2
    * hybridization numbers, m and n
    * hybridizaiton sp^x

    if the number of neighbors > 3, the function only computes the angular 
    defect, the average distances with neighbors and the distance from the 
    average plane.

    Args:
        species (list): list of elements as string
        coords (nat x 3 arry): cartesian coordinates
        rcut (float): cutoff radius to select neighbors

    """

    # compute the distance matrix
    if not distances:
        # compute all distances
        distances = (coords[:, None, :] - coords[None, :, :]) ** 2
        distances = np.sqrt(np.sum(distances, axis=-1))

    # dict of data to set up the dataframe
    natom = len(species)
    columns = ["atom index", "species", "x", "y", "z", "angular defect",
               "pyrA", "improper", "spherical curvature", 
               "c_pi^2", "lambda_pi^2", "m", "n", "hybridization",
               "dist. from ave. plane", "n_neighbors", "ave. neighb. dist."]
    data = pd.DataFrame(columns=columns)
    # data = {
    #     "atom index": np.arange(1, natom + 1).tolist(),
    #     "species": species,
    #     "x": coords[:, 0],
    #     "y": coords[:, 1],
    #     "z": coords[:, 2]
    # }
    # data.update({c: np.full(natom, np.nan) for c in columns[5:]})

    for iat in range(natom):
        # atom A
        atom_a = coords[iat, :]
        line = {"atom index": iat + 1, "species": species[iat]}
        line.update({x: coords[iat, i] for i, x in enumerate("xyz")})

        # look for atoms of *(A)
        n_neighbors = 0
        ave_distance = 0
        star_a = list()
        for jat in range(natom):
            if jat == iat:  # skip current atom
                continue

            rc = get_bond_cutoff(species[iat], species[jat], rcut)
            if distances[iat, jat] < rc:
                star_a.append(coords[jat, :])
                n_neighbors += 1
                ave_distance += distances[iat, jat]
        
        star_a = np.array(star_a, dtype=np.float)

        # number of neighbors
        line["n_neighbors"] = n_neighbors

        # compute data if the number of neighbors is relevant
        if n_neighbors == 3:
            pyrA = curvature.get_pyr_angle(atom_a, star_a)
            line["pyrA"] = pyrA
            line["improper"] = curvature.get_improper(atom_a, star_a)
            line["spherical curvature"] = \
                curvature.get_spherical_curvature(atom_a, star_a)
            
            c_pi, lambda_pi = curvature.get_hybridization_coeff(pyrA)
            line["c_pi^2"] = c_pi ** 2
            line["lambda_pi^2"] = lambda_pi ** 2
            m, n = curvature.get_hybridization_numbers(pyrA)
            line["m"] = m
            line["n"] = n
            line["hybridization"] = curvature.hybridization(pyrA)

        if n_neighbors >= 3:
            line["dist. from ave. plane"] = curvature.get_pyr_distance(atom_a, star_a)
            line["angular defect"] = curvature.get_angular_defect(atom_a, star_a)

        if n_neighbors >= 1:
            line["ave. neighb. dist."] = ave_distance / n_neighbors

        data.loc[iat] = line
    # data["n_neighbors"] = data["n_neighbors"].astype(np.integer, copy=False)
    # df = pd.DataFrame(data, columns=columns)
    # df.set_index(df["atom index"], inplace=True)

    return data, distances


def get_molecule_data(mol, rcut=2.0):
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

    # look for max cutoff
    elements = mol.composition.elements
    bonds = [(elements[iat].symbol, elements[jat].symbol) 
                for iat in range(len(elements)) 
                    for jat in range(iat, len(elements))]
    rcut_max = max([get_bond_cutoff(s1, s2) for s1, s2 in bonds])

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