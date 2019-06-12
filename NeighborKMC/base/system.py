"""### Defines the SystemBase class.

"""


class SystemBase:
    """#### Defines a system class to perform kMC.  
            
    Method assigns an ASE.Atoms object *atoms*
    to the object and assigns a list of sites *sites*.  
    
    A neighbor list *neighbors* is initialized  
    from *sites*, which is checked for inconsistencies  
    by *verify_nlist()*.
    
    **Parameters**  
    *atoms* (ase.atoms, optional): can be passed to connect an ASE atoms  
    object to the system.
    
    *sites* ([Site]): the list of sites that constitute the system.

    **Returns**  
    A SystemBase instance.
    
    **See Also**  
    The module [sites](sites.html)  
    The module [user_sites](../user_sites.html)

    """

    def __init__(self, sites, atoms=None):
        self.sites = sites
        self.neighbors = [s.neighbors for s in sites]
        self.verify_nlist()
        self.atoms = atoms

    def verify_nlist(self):
        """#### Tests the neighborlist for inconsistency.
        
        The method checks if neighborlists are assymetric:  
        If A is a neigbor to B, then is B must  
        also be present in the neighborlist of A.  
        
        **Raises**  
        Warning if the neighborlist is assymetric.
        
        """
        for i, s in enumerate(self.neighbors):
            for nn in s:
                if i not in self.neighbors[nn]:
                    raise Warning("Site " + str(i) + " is a neighbor to site " +
                                  str(nn) + " but not vice-versa")

    def get_ncovs(self, i_site):
        """#### Gets the coverage on nearest neighbor sites.  
            
        Retrieves and returns the occupations of the nearest neighbor
        sites to the site with index *i_site* in *self.sites*.

        **Parameters**  
        *i_site* (int): index of the site relative to *self.sites*.  

        **Returns**  
        *covs* ([int]): list of species occupying the nearest neighbor sites.

            
        """
        covs = [self.sites[n].covered for n in self.neighbors[i_site]]
        return covs

    def find_nn_recurse(self, sim, update_sites, recursion=0):
        """#### Deep serach of first nearest neighbors.

        Calculates the first nearest neighbors to *update\_sites*.

        For example, when passing *update\_sites*gh = [0,1,2]
        the method returns [0,1,2,NN0 of 0, NN1 of 0, NN0 of 1 ...]

        The method is calling itself recursively until the lattice
        is updated, c.f. the locality of nearest neighbor interactions.

        **Parameters**
        *update_sites* ([int]): the site indices to return neighborlist of.

        *recursion* (int, optional): the recursive level of
        which function was called.

        """
        out = [n for n in update_sites]

        for s in update_sites:
            out.extend(self.neighbors[s])

        out = list(set(out))

        if recursion < sim.nninter - 1:
            out = self.find_nn_recurse(sim, out, recursion + 1)

        return out
