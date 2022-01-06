import base.sites


# constant site types
SITE_FCC111_TOP = 0
SITE_FCC111_BRIDGE = 1
SITE_FCC111_FCC_HOLLOW = 2
SITE_FCC111_HCP_HOLLOW = 3


# Species constants
SPECIES_X = 0
SPECIES_COX = 1
SPECIES_OX = 2


class Site(base.sites.SiteBase):
    def __init__(
        self,
        stype=SITE_FCC111_TOP,
        covered=0,      # covered is the species index
        ind=[],         # ase indices of site-related atoms
        lattice_pos=None
    ):
        base.sites.SiteBase.__init__(
            self,
            stype=stype,
            covered=covered,
            ind=ind,
            lattice_pos=lattice_pos
        )


# class TopSite(base.sites.SiteBase):
#     def __init__(
#         self,
#         stype=SITE_FCC111_TOP,
#         covered=0,      # covered is the species index
#         ind=[],         # ase indices of site-related atoms
#         lattice_pos=None
#     ):
#         base.sites.SiteBase.__init__(
#             self,
#             stype=stype,
#             covered=covered,
#             ind=ind,
#             lattice_pos=lattice_pos
#         )


# class BridgeSite(base.sites.SiteBase):
#     def __init__(
#         self,
#         stype=SITE_FCC111_BRIDGE,
#         covered=0,
#         ind=[],         # ase indices of site-related atoms
#         lattice_pos=None
#     ):
#         base.sites.SiteBase.__init__(
#             self,
#             stype=stype,
#             covered=covered,
#             ind=ind,
#             lattice_pos=lattice_pos
#         )


# class FCCHollowSite(base.sites.SiteBase):
#     def __init__(
#         self,
#         stype=SITE_FCC111_FCC_HOLLOW,
#         covered=0,
#         ind=[],         # ase indices of site-related atoms
#         lattice_pos=None
#     ):
#         base.sites.SiteBase.__init__(
#             self,
#             stype=stype,
#             covered=covered,
#             ind=ind,
#             lattice_pos=lattice_pos
#         )


# class HCPHollowSite(base.sites.SiteBase):
#     def __init__(
#         self,
#         stype=SITE_FCC111_HCP_HOLLOW,
#         covered=0,
#         ind=[],         # ase indices of site-related atoms
#         lattice_pos=None
#     ):
#         base.sites.SiteBase.__init__(
#             self,
#             stype=stype,
#             covered=covered,
#             ind=ind,
#             lattice_pos=lattice_pos
#         )
