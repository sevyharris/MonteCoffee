"""Script that runs the tutorial for A adsorption and desorption.
 
"""
import numpy as np
from ase.build import fcc100
from .user_sites import Site
from .user_system import System
from .user_kmc import NeighborKMC
from .user_events import (AAdsEvent, ADesEvent)

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
    tend = 10.  # End time of simulation (s)
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
    atoms = fcc100("Pt", a = a, size = (100,10,1) )
    sites = []

    # Create a site for each surface-atom:
    for i in range(len(atoms)):
        sites.append(Site(stype=0,
                          covered=0, ind=[i]))

    # Instantiate a system, events, and simulation.
    # ---------------------------------------------
    p = System(atoms=atoms, sites=sites)

    # Set the global neighborlist based on distances:
    p.set_neighbors(Ncutoff)

    events = [AAdsEvent, ADesEvent]

    # Specify what events are each others' reverse.
    reverse_events = {0: 1}

#    parameters = {"pCO": pCO, "pO2": pO2, "T": T,
    parameters = { "Name": "A ads/des Simulation",
                  "reverses ": reverse_events}

    # Instantiate simulator object.
    sim = NeighborKMC(system=p, tend=tend,
                      parameters=parameters,
                      events=events,
                      rev_events=reverse_events)

    # Run the simulation.
    sim.run_kmc()
    print("Simulation end time reached ! ! !")

    # Plot the coverages:
    # ------------------------------------------
#    import matplotlib.pyplot as plt
#    import matplotlib as mpl
##
#    fs = 20
#    mpl.rcParams['axes.linewidth'] = 3.5
#    mpl.rcParams['mathtext.default'] = 'rm'
#    mpl.rcParams['mathtext.rm'] = 'Arial'
#    mpl.rcParams['font.family'] = 'Arial'
#    mpl.rc("font", size=fs)
##
#    time = np.loadtxt("time.txt")
#    covs = np.loadtxt("coverages.txt")
##
#    cov_A = [sum([1 for val in covs[i] if val == 1]) / float(len(covs[0])) for i in range(len(covs))]
#
#    plt.plot(time, cov_A, 'k-', lw=2, label="NKMC")

#    plt.legend(loc=0, frameon=False, fontsize="small")
#
#    plt.xlabel("time (s)", fontsize=20)
#    plt.ylabel("Coverage", fontsize=20)
#    plt.gca().tick_params(axis='both', labelsize=fs, length=14, width=3.5, which='major', pad=10)
##
#    plt.savefig('cov_out.pdf',bbox_inches='tight')

if __name__ == '__main__':
    run_test()

