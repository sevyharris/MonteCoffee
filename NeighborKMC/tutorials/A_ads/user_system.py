"""Defines the System Class, derived from base.system module.

The System is supposed to be a singleton that
is passed to a singleton NeighborKMC object.

.. seealso::

   Module :py:mod:`NeighborKMC.base.system`
   Module :py:mod:`NeighborKMC.tutorials.A_ads.user_sites`

"""

import numpy as np
from base.system import SystemBase


class System(SystemBase):
    """Class defines a collection of sites and connected atoms.
            
    Calls the base class system.py constructor, 
    sets the global neighborlist from the individual site's
    neighborlist.

    Attributes
    -----------
    atoms: ase.Atoms
        Can (optionally) be passed to connect an ASE.Atoms
        object to the system. This can be useful for visualization
        of the simulation, for example by setting the ase.Atoms tag
        according to coverages.
    sites: list(Site)
        The sites that constitute the system.

    .. seealso::

       Module :py:mod:`NeighborKMC.base.system`

    """

    def __init__(self, atoms=None, sites=[]):
        SystemBase.__init__(self, atoms=atoms, sites=sites)

    def set_neighbors(self, Ncutoff, pbc=False):
        """Sets neighborlists of self.sites by distances.

        Loops through the sites and using self.atoms, the
        method adds neighbors to the sites that are within a
        neighbor-distance (Ncutoff).

        Parameters
        -----------
        Ncutoff: float
            The cutoff distance for nearest-neighbors in angstroms
        pbc: bool
            If the neighborlist should be computed with periodic boundary
            conditions. To make a direction aperiodic, introduce a vacuum larger
            than Ncutoff in this direction in self.atoms.

        Raises
        ---------
        Warning:
            If self.atoms is not set, because then distances cannot
            be used to determine neighbors.

        """
        if self.atoms is None:
            raise Warning("Tried to set neighbor-distances in user_system.set_neighbors() with self.atom = None")

        # Set the neighbor list for each site using distances.
        # ------------------------------------------
        for i, s in enumerate(self.sites):
            # Position of site
            for j, sother in enumerate(self.sites):
                # Length of distance vector:
                dcur = self.atoms.get_distance(s.ind, sother.ind, mic=pbc)

                # If the site is a neighbor:
                if dcur < Ncutoff and j != i:
                    s.neighbors.append(j)

        if len(self.neighbors) == 0:
            self.neighbors = [s.neighbors for s in self.sites]
            self.verify_nlist()

    def cover_system(self, species, coverage):
        """Covers the system with a certain species.
            
        Randomly covers the system with a species, at a
        certain fractional coverage.
    
        Parameters
        ----------
        species: int
            The species as defined by the user (e.g. empty=0,A=1).
        coverage: float
            The fractional coverage to load lattice with.

        """
        n_covered = int(np.round(coverage * len(self.system.sites)))
        chosen_sites = np.random.choice(len(self.system.sites), n_covered)
        for c in chosen_sites:
            self.system.sites[c].covered = species
