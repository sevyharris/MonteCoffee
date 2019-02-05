# Script that tests MC Code
from __future__ import print_function, division
import sys
import numpy as np
from ase import Atom,Atoms
from ase.visualize import view
from ase.cluster import Octahedron
from user_sites import Site
from user_system import System
from user_kmc import NeighborKMC
from user_constants import mCO,mO2,s0CO,s0O,Asite

T=800. # Temperature
pCO = 2E3 # CO pressure
pO2 = 1E3 # O2 pressure
a = 4.00 # Lattice Parameter (only for ase.atoms)


# Clear up old output files:
np.savetxt("time.txt",[])
np.savetxt("coverages.txt",[])
np.savetxt("stype_ev.txt",[])
np.savetxt("stype_ev_other.txt",[])

# Define the system
atoms = Octahedron("Pt",length=12,cutoff=3,latticeconstant=a)
sites = []

# Define a site for each atom that is free with no pre-defined neighbors.

# First find CN of each atom:
CNS = np.zeros(len(atoms))
for i, at in enumerate(atoms):
    pcur = at.position
    dp = np.sqrt([(p[0]-pcur[0])**2.+(p[1]-pcur[1])**2.+
                 (p[2]-pcur[2])**2. for p in atoms.positions])

    CNS[i] = len([val for val in dp if val < a/np.sqrt(2)+0.01\
                  and val >0.0])

# Define surface atoms as non-bulk:
surface_atom_ids = [i for i in range(len(CNS)) if CNS[i]<12]

stypes = {}

for i,k in enumerate(sorted(list(set(CNS)))):
    stypes[k] = i

print("SITE TYPES, ", sorted(list(set(CNS) )))
# Create a site for each surface-atom:
for i,indic in enumerate(surface_atom_ids):
    sites.append(Site(stype=stypes[CNS[indic]],
                 covered=0,ind = [indic]))

# Instantiate a particle
p = System(atoms=atoms,sites=sites)

# Tagged the ase.atoms by site-types:

for s in surface_atom_ids:
    atoms[s].tag = sorted(list(set(CNS))).index(CNS[s])

maxtag = max([a.tag for a in atoms])
for i,b in enumerate(atoms):
    if i not in surface_atom_ids:
        b.tag = maxtag+1


parameters = {"pCO": pCO,"pO2": pO2,"T": T, "Name": "COOx Simulation"}

sim = NeighborKMC(system = p,tend = 1E-8, parameters=parameters)

result = sim.run_kmc()

# Plot the result:
# ----------------------------------
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.mlab as mlab
import matplotlib as mpl
fs = 20 # Font size
mpl.rcParams['axes.linewidth'] = 3.5 #set the value globally
mpl.rcParams['mathtext.default']='rm'
mpl.rcParams['mathtext.rm'] = 'Arial'
mpl.rcParams['font.family'] = 'Arial'
mpl.rc("font",size=fs)

time = np.loadtxt("time.txt")
covs = np.loadtxt("coverages.txt")

cov_CO = [sum([1 for val in covs[i] if val==1])/float(len(covs[0])) for i in range(len(covs))]
cov_O = [sum([1 for val in covs[i] if val==2])/float(len(covs[0])) for i in range(len(covs))]
cov_free = [sum([1 for val in covs[i] if val==0])/float(len(covs[0])) for i in range(len(covs))]

plt.plot(time,cov_CO,'k-',lw=2,label="CO")
plt.plot(time,cov_O,'r-',lw=2,label="O")
plt.plot(time,cov_free,'m--',lw=2,label="free")
plt.legend(loc=0,frameon=False,fontsize="small")
plt.xlabel("time (s)",fontsize=20)
plt.ylabel("Coverage",fontsize=20)
plt.gca().tick_params(axis='both', labelsize=fs,length=14, width=3.5, which='major',pad=10)

plt.show()
# Print out the TOF:
print("Simualtion end time reached")

