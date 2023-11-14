============
Introduction
============

.. contents:: Table of contents

--------------------

pychemcurv is a python package for structural analyzes of molecular systems or 
solid state materials focusing on the local curvature at an atomic scale. The 
local curvature is then used to compute the hybridization of molecular orbitals.

The main features of the library are available from a 
`Plotly/Dash <https://plot.ly/dash/>`_ web application available
here: `pychemcurv.onrender.com/ <https://pychemcurv.onrender.com>`_.
The web-app allows to upload simple xyz files and compute the local geometrical
properties and the hybridization properties. The application source code is available
in a separate repository at `pychemcurv-app <https://github.com/gVallverdu/pychemcurv-app>`_.

.. figure:: img/screenshot.png
    :align: center
    :width: 300

    Pyramidalization angle of a :math:`C_{28}` fullerene mapped on the structure
    with a colorscale.

Features
========

Pychemcurv is divided in two parts. The first one is a standard python package 
which provides two main classes to compute the local curvature at the atomic 
scale and the hybridization of a given atom. Second, a `Plotly/Dash <https://plot.ly/dash/>`_ web 
application is provided in order to perform a geometrical and electronic
analyzes on molecules or materials.

The web application is available at
`pychemcurv.onrender.com/ <https://pychemcurv.onrender.com>`_.
The web-app allows to upload simple xyz files and compute the local geometrical
properties and the hybridization properties. The application source code is available
in a separate repository at `pychemcurv-app <https://github.com/gVallverdu/pychemcurv-app>`_.

Some jupyter notebooks are provided in the ``notebooks/`` folder and present use cases 
of the classes implemented in this package. You can access to these notebooks
online with `binder <https://mybinder.org/>`_.

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/gVallverdu/pychemcurv.git/2020.6.3
    :alt: binder notebooks

Citing pychemcurv
=================

Please, consider to cite the following paper when using either the `pychemcurv`
library or the web application.

Julia Sabalot-Cuzzubbo, Germain Salvato Vallverdu, Didier Bégué and Jacky Cresson
*Relating the molecular topology and local geometry: Haddon’s pyramidalization angle and the Gaussian curvature*, 
J. Chem. Phys. **152**, 244310 (2020).

.. image:: https://img.shields.io/badge/DOI-doi.org%2F10.1063%2F5.0008368-blue
    :target: https://aip.scitation.org/doi/10.1063/5.0008368
    :alt: DOI


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

    pip install .[app]
    # or
    pip install .[viz]
    # or all extras
    pip install .[app, viz]

    # escape square bracket with zsh
    pip install .\[app, viz\]

If you have installed nglview you have to enable the jupyter extension

::

    jupyter-nbextension enable nglview --py --sys-prefix


The files ``requirements.txt`` and ``environment.yml`` are provided to setup
a full environment with all dependencies. Using pip, in a new environment
you can run the following command to install dependencies

::

    pip install -r requirements.txt

or using ``conda`` you can create the new environment and install all
dependencies in one shot by

::

    conda env create -f environment.yml


The name of the new environment is ``curv``.

Do not forget to enable the jupyter nglview extension (see above).


Install in developper mode
--------------------------

In order to install in developper mode, first create an environment
(using one of the provided file for example) and then install using pip

::

    pip install -e .[app, viz]


If you want to build the documentation you also need to install sphinx.
    

Run the web application
=======================

The web application is available in this separate repository: 
`pychemcurv-app https://github.com/gVallverdu/pychemcurv-app <https://github.com/gVallverdu/pychemcurv-app>`_.
The main aim of the application is to use the pychemcurv 
package and visualize the geometrical or chemical atomic quantities mapped on 
the chemical structure of your system.

The application is available online at this address: 
`pychemcurv.onrender.com/ <https://pychemcurv.onrender.com>`_.

Demo video:

.. raw:: html

    <iframe width="560" height="315" src="https://www.youtube.com/embed/q7UO5Gou-lw" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


In order to run the application locally, you have to install all the dependencies
and in particular ``dash`` and ``dash-bio``. You can do that from the files
``requirements.txt`` or ``environment.yml`` or directly using ``pip``.

Then, clone the github repository on your computer

::

    git clone https://github.com/gVallverdu/pychemcurv-app.git


To run the application, change to the ``pychemcurv-app/`` directory and run the 
``app.py`` file.

::

    [user@computer] (curv) > $ cd pychemcurv-app/
    [user@computer] (curv) > $ python app.py
    Running on http://127.0.0.1:8050/
    Debugger PIN: 065-022-191
    * Serving Flask app "app" (lazy loading)
    * Environment: production
    WARNING: This is a development server. Do not use it in a production deployment.
    Use a production WSGI server instead.
    * Debug mode: on

Open the provided url to use the application.

You can switch off the debug mode by setting ``debug=False`` on the last line of 
the ``app.py`` file.

Common error on local execution
-------------------------------

If the application does not start with an error such as:

::

    socket.gaierror: [Errno 8] nodename nor servname provided, or not known


Go to the last lines of the file app.py and comment/uncomment the last
lines to get something that reads

.. code-block:: python

    if __name__ == '__main__':
        app.run_server(debug=True, host='127.0.0.1')
        # app.run_server(debug=False)



Licence and contact
===================

This software was developped at the `Université de Pau et des Pays de l'Adour
(UPPA) <http://www.univ-pau.fr>`_ in the `Institut des Sciences Analytiques et
de Physico-Chimie pour l'Environement et les Matériaux (IPREM)
<http://iprem.univ-pau.fr/>`_ and the `Institut Pluridisciplinaire de Recherches
Appliquées (IPRA) <http://ipra.univ-pau.fr/>`_ and is distributed under the 
`MIT licence <https://opensource.org/licenses/MIT>`_.

**Authors**

* Germain Salvato Vallverdu: `germain.vallverdu@univ-pau.fr <germain.vallverdu@univ-pau.fr>`_
* Julia Sabalot-cuzzubbo `julia.sabalot@univ-pau.fr  <sabalot.julia@univ-pau.fr>`_
* Didier Bégué: `didier.begue@univ-pau.fr <didier.begue@univ-pau.fr>`_
* Jacky Cresson: `jacky.cresson@univ-pau.fr <jacky.cresson@univ-pau.fr>`_


|UPPA| |CNRS| |IPREM|

.. |UPPA| image:: https://www.univ-pau.fr/skins/uppa_cms-orange/resources/img/logoUPPA.png
  :target: https://www.univ-pau.fr/en/home.html
  :height: 75

.. |IPREM| image:: https://annuaire.helioparc.fr/img/2019/11/logo-9.png
  :target: https://iprem.univ-pau.fr/fr/index.html
  :height: 75

.. |CNRS| image:: http://www.cnrs.fr/themes/custom/cnrs/logo.svg
  :target: http://www.cnrs.fr/
  :height: 75

