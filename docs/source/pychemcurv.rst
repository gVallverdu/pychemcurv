============
Core classes
============

.. automodule:: pychemcurv.core

Vertex classes
--------------

VertexAtom class
~~~~~~~~~~~~~~~~

.. autoclass:: pychemcurv.core.VertexAtom
    :members:

TrivalentVertex class
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: pychemcurv.core.TrivalentVertex
    :members:

POAV: Pi-Orbital Axis Vector
----------------------------

POAV stands for :math:`\pi`-Orbital Axis Vector. The definition of this
vector has its origin in the works of R.C. Haddon. The definitions and the
relation between POAV and the local curvature of a molecule using new 
geometrical object such as the angular defect have been established in our
recent work [JCP2020]_.

Hereafter, the two classes POAV1 and POAV2 aim to compute quantities related
to the two definitions of the POAV vector.

POAV1
~~~~~

.. autoclass:: pychemcurv.core.POAV1
    :members:

POAV2
~~~~~

.. autoclass:: pychemcurv.core.POAV2
    :members:

