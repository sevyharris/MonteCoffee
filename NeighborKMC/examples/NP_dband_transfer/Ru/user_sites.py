"""Defines the Site Class derived from base.site.SiteBase.

The site class is defined here as an interface to the base
class in base.site.SiteBase, where the user can add custom tags.
Custom tags can be used to evaluate the rate-constants in user_events.py

See Also
---------
Module: user_events

"""
from base.sites import SiteBase


class Site(SiteBase):
    """A site object.
           
    Method calls the base class constructor first.  
    Then the user can attach custom variables to site  
    objects, such as coordination numbers, positions, etc.
    
    Attributes
    -------------
    stype: int
        The site type, user must decide what that implies.
        Example: 0 ~ (111) facet ontop, 1 ~ Edge ontop ...

    covered: int
        The species that covers the site, user must decide what the integer implies.
        Example: 0 ~ empty-site, 1 = Oxygen covered, 2 ~ CO covered.

    ind: list(int)
        The atomic-indices c.f. an ASE.Atoms object that constitute
        the site. This is can be used later for visualization.

    See Also
    -----------
    Module: base.sites

    """

    def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None, GCN=12, strain=0.0):
        SiteBase.__init__(self, stype=stype, covered=covered, ind=ind,
                          lattice_pos=lattice_pos)
        self.GCN=GCN
        self.strain=strain
