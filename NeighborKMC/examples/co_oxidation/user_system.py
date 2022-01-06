import numpy as np
import base.system


class System(base.system.SystemBase):
    """Defines the collection of sites and neighbors"""
    def __init__(self, atoms=None, sites=[]):
        base.system.SystemBase.__init__(self, atoms=atoms, sites=sites)

    def set_neighbors(self, Ncutoff, pbc=False):
        if self.atoms is None:
            raise Warning("No ase atoms defined")

        for i, site in enumerate(self.sites):
            for j, other_site in enumerate(self.sites):
                dist = self.atoms.get_distance(site.ind, other_site.ind, mic=pbc)

                if dist < Ncutoff + 0.00001 and j != i:
                    site.neighbors.append(j)

    def cover_system(self, species, coverage):
        # TODO multispecies coverage
        n_covered = int(np.round(coverage * len(self.system.sites)))
        chosen_sites = np.random.choice(len(self.system.sites), n_covered)
        for c in chosen_sites:
            self.system.sites[c].covered = species
