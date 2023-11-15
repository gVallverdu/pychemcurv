# coding: utf-8

"""
The ``pychemcurv.vis`` module implements the ``CurvatureViewer`` 
class in order to visualize a molecule or a periodic structure in a jupyter 
notebook and map a given properties on the atoms using a color scale.

This class needs, `nglview <https://github.com/arose/nglview>`_ and uses
ipywidgets in a jupyter notebook to display the visualization. Run the 
following instructions to install nglview and achieve the configuration 
in order to be able to use nglview in a jupyter notebook

::

    conda install nglview -c conda-forge
    jupyter-nbextension enable nglview --py --sys-prefix

or

::

    pip install nglview
    jupyter-nbextension enable nglview --py --sys-prefix

"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from pymatgen.core import Molecule, Structure
from .analysis import CurvatureAnalyzer

__all__ = ["CurvatureViewer"]


class CurvatureViewer:
    """ This class provides a constructor for a NGLView widget in order to
    visualize the wanted properties using a color scale mapped on the 3D structure
    of the molecule or the structure.
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

        # compute data from CurvatureAnalyzer
        self.data = CurvatureAnalyzer(
            structure, bond_tol, rcut, bond_order).data

    def get_view(self, representation="ball+stick", radius=0.25, aspect_ratio=2,
                 unitcell=False, width="700px", height="500px"):
        """ Set up a simple NGLView widget with the ball and stick or
        licorice representation of the structure.

        Args:
            representation (str): representation: 'ball+stick' or 'licorice'
            radius (float): bond (stick) radius
            aspect_ratio (float): ratio between the balls and stick radiuses
            unitcell (bool): If True and structure is periodic, show the unitcell.
            width (str): width of the nglview widget, default '700px'
            height (str): height of the nglview widget, default '500px'

        Returns:
            Return a ``NGLWidget`` object
        """

        # try to import nglview
        try:
            import nglview as nv
        except ImportError as e:
            print("WARNING: You need to install ase and nglview to perform "
                  "visualization.")
            print(e)
            return None

        if representation not in ["ball+stick", "licorice"]:
            print("Switch representation to 'ball+stick'")

        view = nv.show_pymatgen(self.structure)
        view.clear()
        view.center()
        view.add_representation(
            representation,
            radius=radius,
            aspect_ratio=aspect_ratio,
        )

        # check unitcell
        if isinstance(self.structure, Structure) and unitcell:
            view.add_unitcell()

        # resize nglview widget
        view._remote_call("setSize", target="Widget", args=[width, height])

        return view

    def map_view(self, prop, radius=0.25, aspect_ratio=2, unitcell=False,
                 cm="viridis", minval=None, maxval=None, orientation="vertical",
                 label=None, width="700px", height="500px"):
        """ Map the given properties on a color scale on to the molecule using
        a ball and stick representations. The properties can be either the name
        of a column of the data computed using the CurvatureAnalyzer class, or, 
        an array of values of a custum property. In the last case, the size of
        the array must be consistent with the number of atoms in the system.

        Args:
            prop (str or array): name of the properties or values you want to map
            radius (float): bond (stick) radius
            aspect_ratio (float): ratio between the balls and stick radiuses
            unitcell (bool): If True and structure is periodic, show the unitcell.
            cm (str): colormap from ``matplotlib.cm``.
            minval (float): minimum value to consider for the color sacle
            maxval (float): maximum value to consider for the color sacle
            orientation (str): orientation of the colorbar ``'horizontal'`` or ``'vertical'``
            label (str): Name of the colorbar. If None, use prop.
            width (str): width of the nglview widget, default '700px'
            height (str): height of the nglview widget, default '500px'

        Returns:
            Returns an ipywidgets ``HBox`` or ``VBox`` with the ``NGLWidget``
            and a color bar associated to the mapped properties. The 
            ``NGLWidget`` is the first element of the children, the colorbar
            is the second one.
        """

        # try to import ipywidgets
        try:
            from ipywidgets import HBox, VBox, Output
        except ImportError as e:
            print("You need ipywidgets available with jupyter notebook.")
            print(e)
            return None

        # check property data
        if isinstance(prop, str):
            if prop in self.data.columns:
                prop_vals = self.data[prop].values
                label = prop if label is None else label
            else:
                print("Available data are", data.columns)
                raise ValueError("prop %s not found in data." % prop)
        else:
            try:
                prop_vals = np.array(prop, dtype=np.float64).reshape(
                    len(self.structure))
            except ValueError:
                print("property = ", prop)
                raise ValueError(
                    "Cannot convert prop in a numpy array of floats.")

            # colorbar label
            label = "" if label is None else label

        # check orientation
        if orientation not in ["vertical", "horizontal"]:
            orientation = "horizontal"

        # find property boundary
        if minval is None:
            minval = np.nanmin(prop_vals)
        if maxval is None:
            maxval = np.nanmax(prop_vals)

        # normalize colors
        normalize = mpl.colors.Normalize(minval, maxval)
        cmap = mpl.cm.get_cmap(cm)
        # set up a matplotlib figure for the colorbar
        if orientation == "horizontal":
            _, ax = plt.subplots(figsize=(8, 1))
        else:
            _, ax = plt.subplots(figsize=(1, 8))

        mpl.colorbar.ColorbarBase(ax, cmap=cmap, norm=normalize,
                                  orientation=orientation)
        ax.set_title(label)

        # set up the visualization
        view = self.get_view(representation="ball+stick", radius=radius,
                             aspect_ratio=aspect_ratio, unitcell=unitcell,
                             width=width, height=height)

        # resize nglview widget
        view._remote_call("setSize", target="Widget", args=[width, height])

        # set the atom colors
        for iat, val in enumerate(prop_vals):
            if np.isnan(val):
                continue
            color = mpl.colors.rgb2hex(cmap(X=normalize(val), alpha=1))
            view.add_representation('ball+stick', selection=[iat], color=color,
                                    radius=1.05 * radius,
                                    aspect_ratio=aspect_ratio)

        # resize nglview widget
        view._remote_call("setSize", target="Widget", args=[width, height])

        # place the colorbar in an Output() widget
        out = Output()
        with out:
            plt.show()

        # gather the view and colorbar in a vbox or hbox depending on
        # the orientation
        if orientation == "vertical":
            box = HBox(children=[view, out])
        else:
            box = VBox(children=[view, out])

        return box
