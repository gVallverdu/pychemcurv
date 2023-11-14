# coding: utf-8

"""
This python packages provides classes in order to compute the local curvature
in a molecule or a material at the atomic scale and the hybridization of the
molecular orbitals of the atoms. The `utils` module allows to compute all
the quantities for all the atoms of a molecule or a unit cell.
"""

__author__ = "Germain Salvato-Vallverdu"
__copyright__ = "University of Pau and Pays Adour"
__version__ = "2020.4.22"
__email__ = "germain.vallverdu@univ-pau.fr"

from .core import VertexAtom, TrivalentVertex, POAV1, POAV2
from .analysis import CurvatureAnalyzer
from .vis import CurvatureViewer

# import convenient object from pymatgen
from pymatgen.core import Molecule, Structure
