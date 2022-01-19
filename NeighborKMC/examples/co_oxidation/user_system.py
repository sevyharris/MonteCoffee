import enum
import numpy as np
import base.system


def x_sort(site):
    return site.lattice_pos[0]


def y_sort(site):
    return site.lattice_pos[1]


class System(base.system.SystemBase):
    """Defines the collection of sites and neighbors"""
    def __init__(self, atoms=None, sites=[]):
        base.system.SystemBase.__init__(self, atoms=atoms, sites=sites)

    def set_neighbors(self, Ncutoff, pbc=False):
        if self.atoms is None:
            raise Warning("No ase atoms defined")

        # # this currently scales by n^2, which is terrible
        # for i, site in enumerate(self.sites):
        #     for j, other_site in enumerate(self.sites):
        #         dist = self.atoms.get_distance(site.ind, other_site.ind, mic=pbc)

        #         if dist < Ncutoff + 0.00001 and j != i:
        #             site.neighbors.append(j)

        for i, site in enumerate(self.sites):
            # do a left hand search of the array - index counts down to 0
            for j in range(i - 1, -1, -1):
                other_site = self.sites[j]
                dist = self.atoms.get_distance(site.ind, other_site.ind, mic=pbc)
                if dist < Ncutoff + 0.00001:
                    site.neighbors.append(j)
                if self.atoms[site.ind].position[0] - self.atoms[other_site.ind].position[0] > Ncutoff + 0.00001:
                    # if np.abs(self.atoms[site.ind].position[0] - self.atoms[other_site.ind].position[0]) > 2.0 * Ncutoff:
                    break
            # do a right hand search of the array
            for j in range(i + 1, len(self.sites)):
                other_site = self.sites[j]
                dist = self.atoms.get_distance(site.ind, other_site.ind, mic=pbc)
                if dist < Ncutoff + 0.00001:
                    site.neighbors.append(j)
                if self.atoms[other_site.ind].position[0] - self.atoms[site.ind].position[0] > Ncutoff + 0.00001:
                    # if np.abs(self.atoms[other_site.ind].position[0] - self.atoms[site.ind].position[0]) > 2.0 * Ncutoff:
                    break

    def cover_system(self, species, coverage):
        # TODO multispecies coverage
        n_covered = int(np.round(coverage * len(self.sites)))
        chosen_sites = np.random.choice(len(self.sites), n_covered)
        for c in chosen_sites:
            self.sites[c].covered = species
