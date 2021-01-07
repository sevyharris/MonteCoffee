"""Script that runs a very simple NN-kMC from one file only to have a quick start.

"""

from base.events import EventBase
from ase.build import fcc111
from base.sites import SiteBase
import numpy as np
from base.kmc import NeighborKMCBase
from base.system import SystemBase
from base.logging import Log

class Adsorption(EventBase):

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):

        if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1. * self.params["pA"]
        return R  

    def do_event(self, system, site, other_site):
        # Cover the two sites with species 1
        system.sites[site].covered = 1
        system.sites[other_site].covered = 1

    def get_involve_other(self):
        return True

class Desorption(EventBase):

    def __init__(self, params):
        EventBase.__init__(self, params)

    def possible(self, system, site, other_site):

        if system.sites[site].covered == 1 and system.sites[other_site].covered == 1:
            return True
        else:
            return False

    def get_rate(self, system, site, other_site):
        R = 1.
        return R  

    def do_event(self, system, site, other_site):
        # empty the sites:
        system.sites[site].covered = 0
        system.sites[other_site].covered = 0

    def get_involve_other(self):
        return True

class simple_NKMC(NeighborKMCBase):

    def __init__(self, system, tend, parameters={}, events=[]):
        self.events = [ev(parameters) for ev in events]
        NeighborKMCBase.__init__(self, system=system,
                                 tend=tend, parameters=parameters)

    def run_kmc(self):
    
        logparams = {}
        logparams.update(self.parameters)
        logparams.update({"tend": self.tend,
                           "Nsites": self.system.Nsites,
                           "Number of events": len(self.events),
                           "Number of site-types (stypes)": len(list(set([m.stype for m in self. system.sites])))
                           })
        log = Log(logparams)
    
        stepN_CNT = 0
        stepNMC = 0
     
        while self.t < self.tend:
            self.frm_step()
      
            if stepN_CNT >= self.LogSteps:
                print("Time : ", self.t, "\t Covs :", self.system.get_coverages(self.Nspecies))
                log.dump_point(stepNMC, self.t, self.evs_exec)
                stepN_CNT = 0
    
            stepN_CNT += 1
            stepNMC += 1 
        


########### START the program definitions ################

events = [Adsorption, Desorption]

## Define sites

a0 = 4.00  # Lattice Parameter (not related to DFT!)
atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
sites = []
# Define a site for each atom that is empty with no pre-defined neighbors:
for i in range(len(atoms)):
    sites.append(SiteBase(stype=0, covered=0, ind=[i]))

## Init system, neighborlists
p = SystemBase(atoms=atoms, sites=sites)
Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff

for i, s in enumerate(sites):
      for j, sother in enumerate(sites):
          dcur = atoms.get_distance(s.ind[0], sother.ind[0], mic=True)
          if dcur < Ncutoff and j != i:
              s.neighbors.append(j)

# Init NeighborKMC object
parameters = {"pA": 10., "Name": "Quickstart simulation"}
sim = simple_NKMC(system=p,
                  tend=10.0, # end after 10.s.
                  parameters=parameters, # parameters for event rate-constants.
                  events=events) # the list of events

if __name__ == '__main__':
    sim.run_kmc()
