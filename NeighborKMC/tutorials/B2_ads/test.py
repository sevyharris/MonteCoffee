"""Script that runs the tutorial of B2 adsorption and desorption on a surface 
 
"""
import numpy as np
from ase.build import fcc100
from .user_sites import Site
from .user_system import System
from .user_kmc import NeighborKMC
from .user_events import (B2AdsEvent, B2DesEvent)

def run_test():
    """Runs the simulation of the dissociativ B2 adsorption and desorption over a surface.

    First, constants are defined and old output files cleared.
    Next, the sites, events, system and simulation objects
    are loaded, and the simulation is performed.

    Last, the results are read in from the generated.txt files and compared with a file containing the results from the mean-field model,
    and plotted using matplotlib.

    """
    # Define constants.
    # ------------------------------------------
    tend = 10.  # End time of simulation (s)
    a = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
    # Clear up old output files.
    # ------------------------------------------
    np.savetxt("time.txt", [])
    np.savetxt("coverages.txt", [])

    # Define the sites from ase.Atoms.
    # ------------------------------------------
    atoms = fcc100("Pt", a = a, size = (10,10,1) )
    sites = []
    atoms.write('surf.traj')
    # Create a site for each surface-atom:
    for i in range(len(atoms)):
        sites.append(Site(stype=0,
                          covered=0, ind=i))

    # Instantiate a system, events, and simulation.
    # ---------------------------------------------
    p = System(atoms=atoms, sites=sites)

    # Set the global neighborlist based on distances:
    p.set_neighbors(Ncutoff,pbc = True)
    
    events = [B2AdsEvent, B2DesEvent]

#    parameters = {"pCO": pCO, "pO2": pO2, "T": T,
    parameters = { "Name": "B2 ads/des Simulation"}

    # Instantiate simulator object.
    sim = NeighborKMC(system=p, tend=tend,
                      parameters=parameters,
                      events=events)

    # Run the simulation.
    sim.run_kmc()
    print("Simulation end time reached ! ! !")

if __name__ == '__main__':
    run_test()

