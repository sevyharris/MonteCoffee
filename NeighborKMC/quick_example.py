
"""Script that runs a minimal working example.

The example simply involves and Adsorption and Desorption reaction
with a constant rate-constant.
"""

import numpy as np
from ase.build import fcc111
from user_sites import Site
from user_system import System
from base.events import EventBase
from user_kmc import NeighborKMC


class Adsorption(EventBase):

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):

        if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1000. * self.params["pA"]
        return self.alpha * R  # alpha important for temporal acceleration.

    def do_event(self, system, site, other_site):
        # Cover it with species 1
        system.sites[site].covered = 1
        system.sites[other_site].covered = 1


class Desorption(EventBase):

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):

        if system.sites[site].covered == 1 and system.sites[other_site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, i_site, other_site):
        R = 101.
        return self.alpha * R  # alpha important for temporal acceleration.

    def do_event(self, system, site, other_site):
        # Cover it with species 1
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0


def run_mve():
    """Runs minimal working example of A adsorption on a nanoparticle.

    """
    # Define constants.
    # ------------------------------------------
    tend = 1.0  # End time of simulation (s)
    a0 = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
    pA = 100.

    # Load events:
    events = [Adsorption, Desorption]
    # Specify what events are eac others' reverse.
    reverse_events = {0: 1}

    # Define the sites from ase.Atoms.
    # ------------------------------------------
    atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
    sites = []

    # **Define a site for each atom that is free with no pre-defined neighbors.**
    # Create a site for each surface-atom:
    for i, indic in enumerate(atoms):
        sites.append(Site(stype=0, covered=0, ind=[i]))

    # Instantiate a system, events, and simulation.
    # ---------------------------------------------
    p = System(atoms=atoms, sites=sites)

    # Set the global neighborlist based on nearest neighbor distances:
    p.set_neighbors(Ncutoff)

    parameters = {"pA": pA,
                  "Name": "Quickstart simulation",
                  "reverses": reverse_events}

    # Instantiate simulator object.
    sim = NeighborKMC(system=p, tend=tend,
                      parameters=parameters,
                      events=events,
                      rev_events=reverse_events)

    # Run the simulation.
    sim.run_kmc()
    print("Simulation end time reached ! ! !")


if __name__ == '__main__':
    run_mve()


