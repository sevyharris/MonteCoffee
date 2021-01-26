"""Script that runs A+B2 reaction over a truncated Octahedron 
 
"""
import numpy as np
from ase.io import write
from ase.cluster import Octahedron
import sys
from .user_sites import Site
from .user_system import System
from .user_kmc import NeighborKMC
from .user_events import (B2AdsEvent, B2DesEvent, 
                         AAdsEvent, ADesEvent, ABreactEvent, ADiffEvent, BDiffEvent)
import os

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
    tend = 40.  # End time of simulation (s)
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
    atoms = Octahedron("Pt", 8, cutoff=3, latticeconstant = a)
    sites = []
    write('trunc_octa.traj', atoms)
    # Create a site for each surface-atom:

    # First find CN of each atom:
    CNS = np.zeros(len(atoms))
    for i, at in enumerate(atoms):
        pcur = at.position
        dp = np.sqrt([(p[0] - pcur[0]) ** 2. + (p[1] - pcur[1]) ** 2. +
                       (p[2] - pcur[2]) ** 2. for p in atoms.positions])
        CNS[i] = len([val for val in dp if 0. < val < a / np.sqrt(2) + 0.01])

    # Define surface atoms as non-bulk.
    surface_atom_ids = [i for i in range(len(CNS)) if CNS[i] < 12]
    for i,indic in enumerate(surface_atom_ids):
        if CNS[indic] == 9:
            sstype = 0
        else:
            sstype = 1
        sites.append(Site(stype=sstype,
                          covered=0, ind=indic))

    # Instantiate a system, events, and simulation.
    # ---------------------------------------------
    p = System(atoms=atoms, sites=sites)

    # Set the global neighborlist based on distances:
    p.set_neighbors(Ncutoff,pbc = True)
    
    events = [B2AdsEvent, B2DesEvent, AAdsEvent,ADesEvent,ABreactEvent,ADiffEvent,BDiffEvent]

#    parameters = {"pCO": pCO, "pO2": pO2, "T": T,
    parameters = { "Name": "AB reaction simulation"}

    # Instantiate simulator object.
    sim = NeighborKMC(system=p, tend=tend,
                      parameters=parameters,
                      events=events)

    # Run the simulation.
    sim.run_kmc()
    print("Simulation end time reached ! ! !")

if __name__ == '__main__':
    run_test()

