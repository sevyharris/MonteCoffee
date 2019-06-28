.. _installation:
.. index:: Installation


Installation
**************

Download
--------------

To download :program:`MonteCoffee` , simply clone the git::

	git clone https://gitlab.com/ChemPhysChalmers/MonteCoffee.git

or download a .zip archive::

    https://gitlab.com/ChemPhysChalmers/MonteCoffee/-/archive/master/MonteCoffee-master.zip

Setup
------

:program:`MonteCoffee` is presently installed, simply by copying all
contents of `NeighborKMC/` to an appropriate directory.

To verify that the code functions on your current workstation, you
can run::

	python3 test.py

Dependencies
--------------

:program:`MonteCoffee` depends on the following `Python` packages:
 - `Numpy <https://www.numpy.org/>`_ used for array manipulations.
 - `Matplotlib <https://matplotlib.org/>`_ for plotting results.
 - `ASE <https://wiki.fysik.dtu.dk/ase>`_ for handling atomic structures.
