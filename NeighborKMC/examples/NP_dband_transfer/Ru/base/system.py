"""Defines the SystemBase class.

The module defines a class used to template the System class
defined in user_system.

See Also
---------
Module: sites
Module: user_sites

"""


class SystemBase:
    """Defines a system class to perform kMC.
            
    Method assigns an ASE.Atoms object (atoms)
    to the object and assigns a list of sites (sites).
    
    A neighbor list (neighbors) is initialized
    from the sites, which is checked for inconsistencies
    by the method verify_nlist().
    
    Attributes
    -----------
    atoms: ase.atoms
        Can be passed to connect an ASE atoms  object to the system.
    sites: Site
        The list of sites that constitute the system.


    """

    def __init__(self, sites, atoms=None):
        self.sites = sites
        self.neighbors = [s.neighbors for s in sites]
        self.verify_nlist()
        self.atoms = atoms

    def verify_nlist(self):
        """Tests the neighborlist for inconsistency.
        
        The method checks if neighborlists are assymetric:  
        If A is a neighbor to B, then B must
        also be present in the neighborlist of A.  
        
        Raises
        ---------
        Warning:
            If the neighborlist is assymetric.
        
        """
        for i, s in enumerate(self.neighbors):
            for nn in s:
                if i not in self.neighbors[nn]:
                    raise Warning("Site " + str(i) + " is a neighbor to site " +
                                  str(nn) + " but not vice-versa")

    def get_ncovs(self, i_site):
        """Gets the coverage on nearest neighbor sites.
            
        Retrieves and returns the occupations of the nearest neighbor
        sites to the site with index `i_site` in `self.sites`.

        Parameters
        -----------
        i_site: int
            Index of the site in `self.sites`.

        Returns
        -----------
        covs: list(int)
            List of species occupying the nearest neighbor sites.
            
        """
        covs = [self.sites[n].covered for n in self.neighbors[i_site]]
        return covs

    def find_nn_recurse(self, sim, update_sites, recursion=0):
        """Deep search of first nearest neighbors.

        Calculates the first nearest neighbors for a list of site_indices (update_sites).

        For example, when passing update_sites = [0,1,2],
        the method returns [0,1,2,N neighbor 0 of site 0, Neighbor 1 of site 0, ...,
        Neighbor 0 of site 1, ...].

        The method is calling itself recursively until the lattice
        is updated, c.f. the locality of nearest neighbor interactions.

        Parameters
        ------------
        update_sites: list(int)
            The site indices to return neighborlist of.
        recursion: int
            The recursive level of which function was called, because the method
            calls itself, for example in base.kmc.frm_update().

        Returns
        --------
        out: list(int)
            An update to the list update_sites where the neighbors to update_sites
            are added.

        See Also
        ---------
        kmc.NeighborKMC.frm_update()

        """
        out = [n for n in update_sites]

        for s in update_sites:
            out.extend(self.neighbors[s])

        out = list(set(out))

        if recursion < sim.nninter - 1:
            out = self.find_nn_recurse(sim, out, recursion + 1)

        return out
