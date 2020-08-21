#! /usr/bin/env python

"""Script that runs a full example of CO oxidation.
 
"""
import numpy as np
#from ase.parallel import MPI4PY
from ase.io import *
import os
import sys
from ase_functions import (get_surface_atoms,get_dos_of_surface_sides)
from user_sites import Site
from user_system import System
from user_kmc import NeighborKMC
from user_events import (COAdsEvent, CODesEvent,
                         OAdsEvent, ODesEvent,
                         CODiffEvent, ODiffEvent,
                         COOxEvent)


#os.system('rm -r run*')
#world = MPI4PY()
#rank = world.rank # What number simulation copy am I?
#size = world.size # How many total simulation copies?
#rundir = "run_"+str(rank)  # Name of the dir I create

#os.mkdir(rundir) # Create dir
#os.system("cp kMC_options.cfg "+rundir) # copy kMC_options.cfg here
#os.system("cp Pt_561.traj " +rundir)
#os.system("cp DOSCAR " +rundir)
#os.chdir(rundir)





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
    T = 800.  # Temperature
    pCO = 2E3  # CO pressure
    pO2 = 1E3  # O2 pressure
    tend = 1E6  # End time of simulation (s)
    a = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
    Pt_Pt_bulk_dist=3.91761052866/np.sqrt(2.0)
    mu_bulk=-6.86
    d_band_bulk=-3.302
    f_ratio=1.125
    CN_bulk=12.0
    theta_d=1.8
    # Clear up old output files.
    # ------------------------------------------
    np.savetxt("time.txt", [])
    np.savetxt("coverages.txt", [])
    np.savetxt("sid_ev.txt", [])
    np.savetxt("sid_ev_other.txt", [])

    # Define the sites from ase.Atoms.
    # ------------------------------------------
    atoms=read('Pt_561.traj')
    shell_atoms, N_atoms = get_surface_atoms(atomfile='Pt_561.traj')
    dos_list = get_dos_of_surface_sides(dosfile='DOSCAR',surface_sides=shell_atoms,N_atoms=N_atoms)
    sites = []
    dos_to_GCN_list=[((dd-d_band_bulk)*(2.0*theta_d*f_ratio)/mu_bulk+1)*CN_bulk for dd in dos_list]
    # **Define a site for each atom that is free with no pre-defined neighbors.**
    print ('DOS list')
    print (dos_list)
    print ('transfer GCN list')
    print (dos_to_GCN_list)
    print (max(dos_to_GCN_list),min(dos_to_GCN_list))
    sys.exit() 
    stypes = {}  # site-types: (111),(100),edge,corner.
    for i, k in enumerate(list(set(dos_to_GCN_list))):
        stypes[k] = i

    # Create a site for each surface-atom:
    for i, indic in enumerate(shell_atoms):
        sites.append(Site(stype=stypes[dos_to_GCN_list[i]],
                          covered=0, ind=[indic], GCN=dos_to_GCN_list[i], strain=0.0))

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

#    # Plot the coverages:
#    # ------------------------------------------
#    import matplotlib.pyplot as plt
#    import matplotlib as mpl
#
#    fs = 20
#    mpl.rcParams['axes.linewidth'] = 3.5
#    mpl.rcParams['mathtext.default'] = 'rm'
#    mpl.rcParams['mathtext.rm'] = 'Arial'
#    mpl.rcParams['font.family'] = 'Arial'
#    mpl.rc("font", size=fs)
#
#    time = np.loadtxt("time.txt")
#    covs = np.loadtxt("coverages.txt")
#
#    cov_CO = [sum([1 for val in covs[i] if val == 1]) / float(len(covs[0])) for i in range(len(covs))]
#    cov_O = [sum([1 for val in covs[i] if val == 2]) / float(len(covs[0])) for i in range(len(covs))]
#    cov_free = [sum([1 for val in covs[i] if val == 0]) / float(len(covs[0])) for i in range(len(covs))]
#
#    plt.plot(time, cov_CO, 'k-', lw=2, label="CO")
#    plt.plot(time, cov_O, 'r-', lw=2, label="O")
#    plt.plot(time, cov_free, 'm--', lw=2, label="free")
#    plt.legend(loc=0, frameon=False, fontsize="small")
#    plt.xlabel("time (s)", fontsize=20)
#    plt.ylabel("Coverage", fontsize=20)
#    plt.gca().tick_params(axis='both', labelsize=fs, length=14, width=3.5, which='major', pad=10)
#
#    plt.show()
#
if __name__ == '__main__':
    run_test()


