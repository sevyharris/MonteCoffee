"""Defines the SiteBase Class.

This class is used as a parent for the Site class
defined in user_sites.py.

See Also
---------
Module: user_sites

"""


class SiteBase:
    """Class that templates site objects.

    Assigns a site type (stype) to the site, the
    species covers the site (covered), atomic indices
    that constitute the site (ind), and the sites
    that are nearest-neighbors (neighbors).

    Attributes
    ------------
    stype: int
        The site type. The user must decide what that implies.
        Example: 0 ~ (111)-facet-ontop, 1 ~ edge-ontop ...
    covered: int
        The species that covers the site. The user must decide
        what the integer implies.
        Example: 0 ~ empty-site, 1 ~ Oxygen covered, 2 ~ CO covered.
    ind: list(int)
        The atomic-indices c.f. an ase.Atoms object
        that constitute the site. This is convenient
        to define for later visualization purposes.
    lattice_pos: list(int)
        The lattice position of the site.
        Can be used for systems that obey periodic boundary
        conditions, and to determine neighbor-lists.

    """

    def __init__(self, stype=0, covered=0, ind=[], lattice_pos=None):
        self.stype = stype
        self.covered = covered
        self.ind = ind
        self.lattice_pos = lattice_pos
        self.neighbors = []  # Instantiate with empty neighbor-list
