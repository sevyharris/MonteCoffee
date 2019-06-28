"""Template-script for running a parallel kMC simulation.

The script determines which process number the current machine
(world.rank) is in the set of total processes (world.size).
One directory is created for each process and run separately in
this.

Note
------
This script is to be used in combination with submit_script.sh

"""

import os
import sys
from ase.build import fcc111
from ase.parallel import MPI4PY
from user_sites import Site
from user_system import System
from user_events import A, B, Z
from user_kmc import NeighborKMC
from user_constants import *

world = MPI4PY()
rank = world.rank # What number simulation copy am I?
size = world.size # How many total simulation copies?
rundir = "run_"+str(rank)  # Name of the dir I create

os.mkdir(rundir) # Create dir
os.system("cp kMC_options.cfg "+rundir) # copy kMC_options.cfg here
os.chdir(rundir)


T= float(sys.argv[1])  # Temperature
pA = float(sys.argv[2])  # C2H2 pressure

tend = 1.0  # End time of simulation (s)
a0 = 4.00  # Lattice Parameter (not related to DFT!)
Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff

atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
sites = []

for i, indic in enumerate(atoms):
    sites.append(Site(stype=0, covered=0, ind=[i]))


events = [A, B, Z]
reverse_events = {0: 1}

p = System(atoms=atoms, sites=sites)
p.set_neighbors(Ncutoff)

parameters = {"pA": pA,
              "T": T,
              "Name": "Parallel Simulation",
              "reverse events": reverse_events}


sim = NeighborKMC(system=p,
                  tend=tend,
                  parameters=parameters,
                  events=events,
                  rev_events=reverse_events)

result = sim.run_kmc()
