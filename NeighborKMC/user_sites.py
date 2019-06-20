"""### Defines the Site Class derived from base.site.SiteBase.

The site class is defined here as an interface to the base
class in base.site.SiteBase.

"""
import numpy as np
from base.sites import SiteBase

def set_neighbors(sites, atoms, Ncutoff):
    """Returns a global neighborlist for the sites.

    :param sites:
    :param ncutoff:
    :return:
    """
    # Set the neighbor list for each site using distances.
    # ------------------------------------------
    for i, s in enumerate(sites):
        # Position of site
        pcur = atoms[s.ind[0]].position

        for j, sother in enumerate(sites):
            # Position of potential neighbor site:
            pother = atoms[sother.ind[0]].position

            # Length of distance vector:
            dpabs = np.sqrt((pother[0] - pcur[0]) ** 2. +
                            (pother[1] - pcur[1]) ** 2. +
                            (pother[2] - pcur[2]) ** 2.)

            # If the site is a neighbor:
            if dpabs < Ncutoff and j != i:
                s.neighbors.append(j)


class Site(SiteBase):
    """Constructor for site objects.
           
    Method calls the base class constructor first.  
    Then the user can attach custom variables to site  
    objects, such as coordination numbers, positions, etc.
    
    
    **Parameters**  
    *stype* (int): the site type, user must decide what that implies.  
                 Example: 0 ~ (111) facet ontop, 1 ~ Edge ontop ...  

    *covered* (int): the species that covers the site, user must decide 
                   what the integer implies.  
                   Example: 0 ~ empty-site, 1 = Oxygen covered, 2 ~ CO covered.

    *ind* ([int]): the atomic-indices c.f. an ASE.Atoms object that constitute  
                 the site. This is can be used later for visualization.


    **Returns**  
    A Site instance


    **See Also**  
    The module [base.sites](base/sites.html)

    """

    def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None):
        SiteBase.__init__(self, stype=stype, covered=covered, ind=ind,
                          lattice_pos=lattice_pos)
