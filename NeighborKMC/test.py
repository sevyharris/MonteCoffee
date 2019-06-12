"""#### Script that runs a full example of CO oxidation.

The script sets up a simulation of CO oxidation over Pt using
coordination numbers for the reaction energy landscape. First,    
constants are defined and old output files cleared. 
Next, the sites, events, system and simulation objects  
are loaded, and the simulation is performed.
 
"""
import numpy as np
from ase.cluster import Octahedron
from user_sites import Site
from user_system import System
from user_kmc import NeighborKMC
from user_events import (COAdsEvent, CODesEvent,
                         OAdsEvent, ODesEvent,
                         CODiffEvent, ODiffEvent,
                         COOxEvent)

# Define constants.
# ------------------------------------------
T = 800.  # Temperature
pCO = 2E3  # CO pressure
pO2 = 1E3  # O2 pressure
tend = 1E-7  # End time of simulation (s)
a = 4.00  # Lattice Parameter (not related to DFT!)
Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff

# Clear up old output files.
# ------------------------------------------
np.savetxt("time.txt", [])
np.savetxt("coverages.txt", [])
np.savetxt("stype_ev.txt", [])
np.savetxt("stype_ev_other.txt", [])

# Define the sites from ase.Atoms.
# ------------------------------------------
atoms = Octahedron("Pt", length=12, cutoff=3, latticeconstant=a)
sites = []

# **Define a site for each atom that is free with no pre-defined neighbors.**

# First find CN of each atom:
CNS = np.zeros(len(atoms))
for i, at in enumerate(atoms):
    pcur = at.position
    dp = np.sqrt([(p[0] - pcur[0]) ** 2. + (p[1] - pcur[1]) ** 2. +
                  (p[2] - pcur[2]) ** 2. for p in atoms.positions])

    CNS[i] = len([val for val in dp if 0. < val < a / np.sqrt(2) + 0.01])

# Define surface atoms as non-bulk.
surface_atom_ids = [i for i in range(len(CNS)) if CNS[i] < 12]

stypes = {}  # site-types: (111),(100),edge,corner.
for i, k in enumerate(sorted(list(set(CNS)))):
    stypes[k] = i

# Create a site for each surface-atom:
for i, indic in enumerate(surface_atom_ids):
    sites.append(Site(stype=stypes[CNS[indic]],
                      covered=0, ind=[indic]))

# Set the neighbor list for each site using distances.
# ------------------------------------------
positions = atoms.positions

for i, s in enumerate(sites):
    # Position of site
    pcur = atoms[s.ind[0]].position

    for j, sother in enumerate(sites):
        # Position of potential neighbor site:
        pother = atoms[sother.ind[0]].position

        # Length of distance vector:
        dpabs = np.sqrt((pother[0] - pcur[0]) ** 2. +
                        (pother[1] - pcur[1]) ** 2. +
                        (pother[2] - pcur[2]) ** 2.)

        # If the site is a neighbor:
        if dpabs < Ncutoff and j != i:
            s.neighbors.append(j)

# Instantiate a system and simulation.
# ------------------------------------------

events = [COAdsEvent, CODesEvent, OAdsEvent,
          ODesEvent, CODiffEvent,
          ODiffEvent, COOxEvent]

p = System(atoms=atoms, sites=sites)

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
result = sim.run_kmc()
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
