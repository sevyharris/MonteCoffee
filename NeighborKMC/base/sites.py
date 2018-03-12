r"""
Module: sites.py
The particle module defines the SiteBase Class.

"""

class SiteBase:


    def __init__(self, stype = 0,covered=0, ind=[]):
        r"""Constructor for SiteBase objects.
            
            Method assigns an site type 'stype' to the site, the
            species covers the site 'covered', atomic inidces 
            that constitute the site 'ind', and the sites 
            that are nearest neighbors 'neighbors'.
    
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
            SiteBase instance

        """

        assert type(ind) is list
        self.stype = stype
        self.covered = covered 
        self.ind = ind
        self.neighbors = []
