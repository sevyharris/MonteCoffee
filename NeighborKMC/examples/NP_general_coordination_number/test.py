"""Script that runs a full example of CO oxidation.
 
"""
import numpy as np
from ase.io import write
from user_sites import Site, GCN_setup
from user_system import System
from user_kmc import NeighborKMC
import os,sys
from user_events import (COAdsEvent, CODesEvent,
                         OAdsEvent, ODesEvent,
                         CODiffEvent, ODiffEvent,
                         COOxEvent)


def run_test():
    """Runs the test of CO oxidation over Pt.

    Sets up a simulation of CO oxidation over Pt using
    coordination numbers for the reaction energy landscape.

    First, constants are defined and old output files cleared.
    Next, the sites, events, system and simulation objects
    are loaded, and the simulation is performed.

    Last, the results are read in from the generated.txt files,
    and plotted using matplotlib.

    """
    # Define constants.
    # ------------------------------------------
    T = 700.  # Temperature
    pCO = 2E3  # CO pressure
    pO2 = 1E3  # O2 pressure
    tend = 10.0  # End time of simulation (s)
    a = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff

    # Clear up old output files.
    # ------------------------------------------
    np.savetxt("time.txt", [])
    np.savetxt("coverages.txt", [])
    np.savetxt("evs_exec.txt", [])
    np.savetxt("mcstep.txt", [])

    os.system("rm detail_site_event_evol.hdf5")

    # Define the sites from ase.Atoms.
    # ------------------------------------------
    atoms, GCN_dic, surface_atoms = GCN_setup()

    sites = []

    # Create a site for each surface-atom:  stype corresponse to the GCN
    for GCN_key in GCN_dic.keys():
        for i in GCN_dic[GCN_key]:
          sites.append(Site(stype=float(GCN_key),
                            covered=0, ind=[i], lattice_pos = atoms.positions[i]))
        print (GCN_key)
    
    if not len(surface_atoms) == len(sites):
      print ('Not all surface atoms maped in sites')
      sys.exit()

    # Instantiate a system, events, and simulation.
    # ---------------------------------------------
    p = System(atoms=atoms, sites=sites)

    # Set the global neighborlist based on distances:
    p.set_neighbors(Ncutoff)

    events = [COAdsEvent, CODesEvent, OAdsEvent,
              ODesEvent, CODiffEvent,
              ODiffEvent, COOxEvent]

    # Specify what events are each others' reverse.
    reverse_events = {0: 1, 2: 3, 4: 4, 5: 5}

    parameters = {"pCO": pCO, "pO2": pO2, "T": T,
                  "Name": "COOx Simulation",
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


