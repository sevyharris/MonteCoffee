r"""
Module: user_sites.py

The particle module defines the Site Class derived 
from base.site.SiteBase.
"""

from base.sites import SiteBase

class Site(SiteBase):
    r"""Constructor for Site objects.
            
            Method simply calls the base class constructor.
    
            Parameters
            ----------
            stype : int
                The site type, user must decide what that implies.
                Example: 0 ~ (111) facet ontop, 1 ~ Edge ontop ...

            covered : int
                The species that covers the site, user must decide
                what the integer implies.
                Example: 0 ~ empty-site, 1 = Oxygen covered,
                2 ~ CO covered.

            ind : list of int
                The atomic-indices c.f. an ASE.Atoms object 
                that constiture the site. This is convenient
                to define for later visualization purposes.


            Returns
            -------
            Site instance


            See Also
            --------
            The module base.sites

        """

    def __init__(self, stype = 0,covered=0, ind=[],lattice_pos=None):
        SiteBase.__init__(self,stype=stype,covered=covered,ind=ind,
                          lattice_pos=lattice_pos)
