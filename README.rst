==========
pychemcurv
==========

.. image:: https://readthedocs.org/projects/pychemcurv/badge/?version=latest
    :target: https://pychemcurv.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Status

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/gVallverdu/pychemcurv.git/2020.6.3
    :alt: binder notebooks

.. image:: https://img.shields.io/badge/DOI-doi.org%2F10.1063%2F5.0008368-blue
    :target: https://aip.scitation.org/doi/10.1063/5.0008368
    :alt: DOI


* `installation <#installation>`_
* `documentation <https://pychemcurv.readthedocs.io/>`_
* Dash web application `source code <https://github.com/gVallverdu/pychemcurv-app>`_, `live demo <https://pychemcurv.onrender.com>`_
* `Notebooks <https://nbviewer.jupyter.org/github/gVallverdu/pychemcurv/tree/master/notebooks/>`_

pychemcurv is a python package for structural analyzes of molecular systems or
solid state materials focusing on the local curvature at an atomic scale. The
local curvature is then used to compute the hybridization of molecular orbitals.

Features
========

Pychemcurv is divided in two parts. The first one is a standard python package
which provides two main classes to compute the local curvature at the atomic
scale and the hybridization of a given atom. Second, a 
`Plotly/Dash <https://plot.ly/dash/>`_ web
application is provided in order to perform geometrical and electronic
analyzes on molecules or materials. The web application is available at
`pychemcurv.onrender.com/ <https://pychemcurv.onrender.com>`_.
The web-app allows to upload simple xyz files and compute the local geometrical
properties and the hybridization properties. The application is available
in a separate repository at `pychemcurv-app <https://github.com/gVallverdu/pychemcurv-app>`_.

Some jupyter notebooks are provided in the ``notebooks/`` folder and present use cases
of the classes implemented in this package. You can access to these notebooks
online with `binder <https://mybinder.org/v2/gh/gVallverdu/pychemcurv.git/2020.6.3>`_.


Citing pychemcurv
=================

Please, consider to cite the following papers when using either the `pychemcurv`
library or the web application.

.. image:: https://img.shields.io/badge/DOI-doi.org%2F10.1063%2F5.0008368-blue
    :target: https://aip.scitation.org/doi/10.1063/5.0008368
    :alt: DOI

Julia Sabalot-Cuzzubbo, Germain Salvato Vallverdu, Didier Bégué and Jacky Cresson
*Relating the molecular topology and local geometry: Haddon’s pyramidalization angle and the Gaussian curvature*, 
J. Chem. Phys. **152**, 244310 (2020).


.. image:: https://img.shields.io/badge/DOI-doi.org%2F10.1063%2F5.0008368-blue
    :target: https://aip.scitation.org/doi/10.1063/5.0170800
    :alt: DOI

Julia Sabalot-Cuzzubbo, N. Cresson, Germain Salvato Vallverdu, Didier Bégué and Jacky Cresson
*Haddon’s POAV2 vs POAV theory for non-planar molecules*, 
J. Chem. Phys. **159**, 174109 (2023).

Installation
============

Before installing pychemcurv it is recommanded to create a virtual environment 
using conda or virtuelenv.

Short installation
------------------

Using pip directly from github, run

::

    pip install git+git://github.com/gVallverdu/pychemcurv.git


Alternatively, you can first clone the pychemcurv repository

:: 

    git clone https://github.com/gVallverdu/pychemcurv.git

and then install the module and its dependencies using

::

    pip install .



Full installation
-----------------

If you want to use the web application locally or if you want to use
`nglview <https://github.com/arose/nglview>`_ to display structures in 
jupyter notebooks you need to install more dependencies. The setup configuration
provides the ``viz`` and ``app`` extras so, using pip, run one of

:: 

    pip install .
    # or
    pip install .[viz]

    # escape square bracket with zsh
    pip install .\[viz\]

If you have installed nglview you have to enable the jupyter extension

::

    jupyter-nbextension enable nglview --py --sys-prefix


The files ``requirements.txt`` and ``environment.yml`` are provided to setup
a full environment with all dependencies. Using pip, in a new environment
you can run

::

    pip install -r requirements.txt

or using ``conda`` you can create the new environment and install all
dependencies in one shot by

::

    conda env create -f environment.yml


The name of the new environment is `curv`.

Do not forget to enable the jupyter nglview extension (see above).


Install in developper mode
--------------------------

In order to install in developper mode, first create an environment
(using one of the provided file for example) and then install using pip

::

    pip install -e .[viz]


If you want to build the documentation you also need to install sphinx.
    

Licence and contact
===================

This software was developped at the `Université de Pau et des Pays de l'Adour
(UPPA) <http://www.univ-pau.fr>`_ in the `Institut des Sciences Analytiques et
de Physico-Chimie pour l'Environement et les Matériaux (IPREM)
<http://iprem.univ-pau.fr/>`_ and the `Institut Pluridisciplinaire de Recherches
Appliquées (IPRA) <http://ipra.univ-pau.fr/>`_ and is distributed under the 
`MIT licence <https://opensource.org/licenses/MIT>`_.


Authors
-------

* Germain Salvato Vallverdu: `germain.vallverdu@univ-pau.fr <germain.vallverdu@univ-pau.fr>`_
* Julia Sabalot-cuzzubbo `julia.sabalot@univ-pau.fr  <sabalot.julia@univ-pau.fr>`_
* Didier Bégué: `didier.begue@univ-pau.fr <didier.begue@univ-pau.fr>`_
* Jacky Cresson: `jacky.cresson@univ-pau.fr <jacky.cresson@univ-pau.fr>`_
