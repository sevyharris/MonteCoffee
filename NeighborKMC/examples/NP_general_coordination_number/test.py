"""Script that runs a full example of CO oxidation.
 
"""
import numpy as np
from ase.io import *
from user_sites import Site
from user_system import System
from user_kmc import NeighborKMC
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
    T = 800.  # Temperature
    pCO = 2E3  # CO pressure
    pO2 = 1E3  # O2 pressure
    tend = 1E-6  # End time of simulation (s)
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
    atoms=read('Pt_561.traj')
    sites = []

    # **Define a site for each atom that is free with no pre-defined neighbors.**
    dist_max=3.0
    pos=atoms.get_positions()
    Natoms=atoms.get_number_of_atoms()
    shell_atoms=[]
    list_GCN=[]
    for ii in xrange(Natoms):
       dists=[atoms.get_distance(int(ii), int(mm)) for mm in xrange(Natoms) if mm not in [ii]]
       dists.sort()
       count = 0
       for kk in dists:
         if kk < dist_max:
           count +=1
         else:
           break
       if count < 12:
         shell_atoms.append([ii,count])
    CN_atoms=[]
    N_shell=len(shell_atoms)
    #print 'N shell atoms', N_shell
    CN_dict={}
    for jj in xrange(N_shell):
      atom_id=shell_atoms[jj][0]
      dists=[]
      for ll in shell_atoms:
        shell_id=ll[0]
        shell_CN=ll[1]
        if not shell_id == atom_id:
          dists.append([atoms.get_distance(atom_id, shell_id),shell_CN])
      dists.sort()
      CN_j=0
      N_surf_shell=0
      for oo in dists:
        if oo[0] < dist_max:
          CN_j+=oo[1]
          N_surf_shell+=1
        else:
          break
      N_bulk=shell_atoms[jj][1]-N_surf_shell
      CN_j += 12.0*N_bulk
    #  print atom_id, shell_atoms[jj][1], N_surf_shell, N_bulk, CN_j, CN_j/12.0
      CN_j /= 12.0
      list_GCN.append(CN_j)

    stypes = {}  # site-types: (111),(100),edge,corner.
    for i, k in enumerate(list(set(list_GCN))):
        stypes[k] = i

    # Create a site for each surface-atom:
    for i, indic in enumerate(shell_atoms):
        sites.append(Site(stype=stypes[list_GCN[indic]],
                          covered=0, ind=[indic]))

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

    # Plot the coverages:
    # ------------------------------------------
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    fs = 20
    mpl.rcParams['axes.linewidth'] = 3.5
    mpl.rcParams['mathtext.default'] = 'rm'
    mpl.rcParams['mathtext.rm'] = 'Arial'
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rc("font", size=fs)

    time = np.loadtxt("time.txt")
    covs = np.loadtxt("coverages.txt")

    cov_CO = [sum([1 for val in covs[i] if val == 1]) / float(len(covs[0])) for i in range(len(covs))]
    cov_O = [sum([1 for val in covs[i] if val == 2]) / float(len(covs[0])) for i in range(len(covs))]
    cov_free = [sum([1 for val in covs[i] if val == 0]) / float(len(covs[0])) for i in range(len(covs))]

    plt.plot(time, cov_CO, 'k-', lw=2, label="CO")
    plt.plot(time, cov_O, 'r-', lw=2, label="O")
    plt.plot(time, cov_free, 'm--', lw=2, label="free")
    plt.legend(loc=0, frameon=False, fontsize="small")
    plt.xlabel("time (s)", fontsize=20)
    plt.ylabel("Coverage", fontsize=20)
    plt.gca().tick_params(axis='both', labelsize=fs, length=14, width=3.5, which='major', pad=10)

    plt.show()

if __name__ == '__main__':
    run_test()


