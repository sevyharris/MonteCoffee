"""### Defines the Site Class derived from base.site.SiteBase.

The site class is defined here as an interface to the base
class in base.site.SiteBase.

"""

from base.sites import SiteBase


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
