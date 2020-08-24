"""Defines the Site Class derived from base.site.SiteBase.

The site class is defined here as an interface to the base
class in base.site.SiteBase, where the user can add custom tags.
Custom tags can be used to evaluate the rate-constants in user_events.py

.. seealso:: Module :py:mod:`NeighborKMC.tutorials.A_ads.user_events`

"""
#See Also
#---------
#Module: user_events

#"""
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
        Here: 0 refers to on top adsorption

    covered: int
        The species that covers the site, user must decide what the integer implies.
        Here: 0 = empty-site and 1 = A covered

    ind: list(int)
        The atomic-indices c.f. an ASE.Atoms object that constitute
        the site. This is can be used later for visualization.


    .. seealso:: Module :py:mod:`NeighborKMC.base.sites`

    """

    def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None):
        SiteBase.__init__(self, stype=stype, covered=covered, ind=ind,
                          lattice_pos=lattice_pos)
