"""### Defines the SiteBase Class.

"""


class SiteBase:
    """#### Class that templates site objects.  
            
    Assigns a site type *stype* to the site, the
    species covers the site *covered*, atomic inidces 
    that constitute the site *ind*, and the sites 
    that are nearest-neighbors *neighbors*.
    
    **Parameters**  
    *stype* (int): the site type. The user must decide what that implies.  
    Example: 0 ~ (111) facet ontop, 1 ~ Edge ontop ...

    *covered* (int): the species that covers the site. The user must decide  
    what the integer implies.  
    Example: 0 ~ empty-site, 1 = Oxygen covered, 2 ~ CO covered.  

    *ind* ([int]): the atomic-indices c.f. an ase.Atoms object   
    that constiture the site. This is convenient  
    to define for later visualization purposes.

    *lattice_pos* ([int]): the lattice position of the site. 
    Can be used for systems that obey periodic boundary
    conditions, and to determine neighbor-lists. 


    **Returns**  
    A SiteBase instance.
    
    **See Also**  
    The module [user_sites](../user_sites.html)

    """

    def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None):
        self.stype = stype
        self.covered = covered
        self.ind = ind
        self.lattice_pos = lattice_pos
        self.neighbors = []  # Instantiate with empty neighbor-list
