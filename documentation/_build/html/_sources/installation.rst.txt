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

:program:`MonteCoffee` is presently installed, simply by copying (or linking) the base-file folder to the
run-directory which contains the user-specific files. 

To verify that the code runs successfully on your current workstation, you
can run::

	python3 test.py

in one of the tutorial folders.

Dependencies
--------------

:program:`MonteCoffee` depends on the following `Python` packages:
 - `Numpy <https://www.numpy.org/>`_ used for array manipulations.
 - `Matplotlib <https://matplotlib.org/>`_ for plotting results.
 - `ASE <https://wiki.fysik.dtu.dk/ase>`_ for handling atomic structures.
 - `h5py <https://docs.h5py.org/en/stable/>`_ for handling file output and compression
 - `mpi4py <https://mpi4py.readthedocs.io/en/stable/#>`_ for running various kMC simulations in parallel 
