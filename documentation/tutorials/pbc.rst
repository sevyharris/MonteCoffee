.. _pbc:
.. index:: Periodic boundary conditions

Periodic boundary conditions
*************************************

Periodic boundary conditions are implemented in :program:`MonteCoffee` in the method `System.set_neighbors() <../api/NeighborKMC.html#NeighborKMC.user_system.System.set_neighbors>`_ if one passes :code:`pbc=True` to the method. The method adds sites to each others neighborlist that
differ in position by a lattice vector plus a nearest neighbor distance: :math:`d_{ij} < |\vec{a}_{\ell}| + d_{N-N}`.

For example, if a Pt(111) surface is instantiated, we can use the `get_distances() <https://wiki.fysik.dtu.dk/ase/_modules/ase/atoms.html#Atoms.get_distances>`_ method from ASE using the `minimum image convention`:

.. code-block:: python

    import numpy as np
    from ase.build import fcc111
    from user_sites import Site
    from user_system import System
    
    
    a0 = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
    
    atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
    sites = []
    
    # Define a site for each atom that is empty with no pre-defined neighbors:
    for i in range(len(atoms)):
        sites.append(Site(stype=0, covered=0, ind=[i]))
        
    p = System(atoms=atoms, sites=sites)

    p.set_neighbors(Ncutoff, pbc=True)
    
 
