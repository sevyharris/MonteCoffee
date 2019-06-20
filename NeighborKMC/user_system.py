"""### Defines the System Class, derived from base.system module.  

*System* is supposed to be a singleton that
is passed to a singleton *NeighborKMC* object.

See Also  
The module [base.system](base/system.html)  
The module [base.kmc](base/kmc.html)  
The module [user_kmc](user_kmc.html)  

"""

import numpy as np
from base.system import SystemBase



class System(SystemBase):
    """#### Class defines a collection of sites and connected atoms.
            
    Calls the base class system.py constructor, 
    sets the global neighborlist from the invididual site
    neighborlist.

    **Parameters**  
    *atoms* (ase.Atoms): Can *optionally* be passed to connect  
    an ASE atoms object to the system.  
    
    *sites* ([Site]): the sites that constitute the system.  


    **Returns**  
    System instance

    **See Also**  
    The module [base.system](base/system.html)  
    The module [base.sites](base/sites.html)  
    The module [user_sites](user_sites.html)  
    
    """

    def __init__(self, atoms=None, sites=[]):
        SystemBase.__init__(self, atoms=atoms, sites=sites)

    def set_neighbors(self, atoms, Ncutoff):
        """Returns a global neighborlist for the sites.

        :param sites:
        :param ncutoff:
        :return:
        """
        # Set the neighbor list for each site using distances.
        # ------------------------------------------
        for i, s in enumerate(self.sites):
            # Position of site
            pcur = atoms[s.ind[0]].position

            for j, sother in enumerate(self.sites):
                # Position of potential neighbor site:
                pother = atoms[sother.ind[0]].position

                # Length of distance vector:
                dpabs = np.sqrt((pother[0] - pcur[0]) ** 2. +
                                (pother[1] - pcur[1]) ** 2. +
                                (pother[2] - pcur[2]) ** 2.)

                # If the site is a neighbor:
                if dpabs < Ncutoff and j != i:
                    s.neighbors.append(j)

    def cover_system(self, species, coverage):
        """#### Covers the system with a certain species.
            
        Randomly covers the system with a species *species*, at a 
        certain fractional coverage *coverage*.
    
        **Parameters**  
        *species* (int): the species as defined by hte user (e.g. empty=0,CO=1).

        *coverage* (float): the fractional coverage to load lattice with.

        """
        n_covered = int(np.round(coverage * len(self.system.sites)))
        chosen_sites = np.random.choice(len(self.system.sites), n_covered)
        for c in chosen_sites:
            self.system.sites[c].covered = species
