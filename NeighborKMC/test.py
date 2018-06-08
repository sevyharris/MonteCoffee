# Script that tests MC Code
import numpy as np
from ase import Atom,Atoms
from ase.visualize import view
from ase.cluster import Octahedron
from user_sites import Site
from user_system import System
from user_kmc import NeighborKMC
from user_constants import mCO,mO2,s0CO,s0O,Asite

T=1100. # Temperature
pCO = 2E3 # CO pressure
pO2 = 1E3 # O2 pressure
a = 4.00 # Lattice Parameter (only for ase.atoms)

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

print "SITE TYPES, ", sorted(list(set(CNS) ))
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


parameters = {"pCO":pCO,"pO2":pO2,"T":T,
              "Name":"COOx Simulation"}

sim = NeighborKMC(system = p,tend = 1E9, parameters=parameters)
result = sim.run_kmc()
