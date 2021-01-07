"""Script that runs a full example of CO oxidation.
 
"""
import numpy as np
from ase.io import write
from ase.build import fcc111
from ase.parallel import MPI4PY
import sys
from user_sites import Site
from user_system import System
from user_kmc import NeighborKMC
from user_events import (OAdsEvent, ODesEvent, 
                         COAdsEvent, CODesEvent, COOxEvent, CODiffEvent, ODiffEvent)
import os

world = MPI4PY()
rank = world.rank # What number simulation copy am I?
size = world.size # How many total simulation copies?
rundir = "run_"+str(rank)  # Name of the dir I create
#os.system('rm -r run_*')
os.mkdir(rundir) # Create dir
os.system("cp kMC_options.cfg "+rundir) # copy kMC_options.cfg here
os.chdir(rundir)


def run_test():
    """Runs the test of A adsorption and desorption over a surface.

    First, constants are defined and old output files cleared.
    Next, the sites, events, system and simulation objects
    are loaded, and the simulation is performed.

    Last, the results are read in from the generated.txt files,
    and plotted using matplotlib.

    """
    # Define constants.
    # ------------------------------------------
    tend = 4.  # End time of simulation (s)
    T = 800.   # Temperature of simulation in K
    pCO = 2E3  # CO pressure in Pa
    pO2 = 1E3  # O2 pressure in Pa
    a = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
    # Clear up old output files.
    # ------------------------------------------
    np.savetxt("time.txt", [])
    np.savetxt("coverages.txt", [])
    np.savetxt("sid_ev.txt", [])
    np.savetxt("sid_ev_other.txt", [])

    # Define the sites from ase.Atoms.
    # ------------------------------------------
    atoms = fcc111("Pt", a = a, size = (10,10,1))
    atoms.write('surface.traj')
    sites = []
    # Create a site for each surface-atom:
    for i in range(len(atoms)):
        sites.append(Site(stype=0,
                          covered=0, ind=[i]))

    # Instantiate a system, events, and simulation.
    # ---------------------------------------------
    p = System(atoms=atoms, sites=sites)

    # Set the global neighborlist based on distances:
    p.set_neighbors(Ncutoff,pbc = True)
    
    events = [COAdsEvent, CODesEvent, OAdsEvent, ODesEvent, CODiffEvent, ODiffEvent, COOxEvent]

    # Specify what events are each others' reverse.
    reverse_events = {0:1,2:3,4:4,5:5}
    #reverse_events = {0:1, 2:3, 4:4, 5:5}

    parameters = {"pCO": pCO, "pO2": pO2, "T": T,
                  "Name": "COOx Pt(111) reaction simulation",
                  "reverses ": reverse_events}

    # Instantiate simulator object.
    sim = NeighborKMC(system=p, tend=tend,
                      parameters=parameters,
                      events=events,
                      rev_events=reverse_events)

    # Run the simulation.
    sim.run_kmc()
    print("Simulation end time reached ! ! !")

if __name__ == '__main__':
    run_test()

