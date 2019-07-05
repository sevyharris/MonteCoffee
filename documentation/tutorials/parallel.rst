.. _parallel:
.. index:: Parallelisation

Parallel simulations
*************************************
kMC simulations are stochastic in nature, making it reasonable to perform multiple identically prepared simulations
to assess the convergence. In the following, we assume that `ASE <https://wiki.fysik.dtu.dk/ase/>`_ and `MPI4PY <https://pypi.org/project/mpi4py/>`_ are installed.

To submit a simulation, defined in a file `kmc_master_parallel.py`, in an environment that implements `SLURM <https://slurm.schedmd.com/>`_, one can submit a bash script as

.. code-block:: bash

    #!/usr/bin/env bash
    #SBATCH -N 1
    #SBATCH -n 20
    #SBATCH -p PROJNAME
    #SBATCH -A ACCNO
    #SBATCH -o out.txt
    #SBATCH -e err.txt
    #SBATCH -t 96:15:00
    #SBATCH -J project_dir/submit_dir
    #SBATCH --mail-user=USER@UNI.se
    #SBATCH --mail-type=END
 
    export NCORES=$(($SLURM_JOB_NUM_NODES*16))

    cp *.py $TMPDIR
    cp kMC_options.cfg $TMPDIR
    cp -r base $TMPDIR

    cd $TMPDIR


    while sleep 1800; do
        # This will be executed once per every 3600 seconds
        cp -r run_*/ $SLURM_SUBMIT_DIR
    done &     # The &-sign after the done-keyword places 
               # the while-loop in a sub-shell in the background
    LOOPPID=$! # Save the PID of the subshell running the loop


    mpirun -np 10 python kmc_master_parallel.py 273. 1E2 1E3 1E3 18 18

    # Copy the files back after run:

    kill $LOOPPID

    cp -r run_*/ $SLURM_SUBMIT_DIR
    
    
This script copies all python files in MonteCoffee to a compute node, and cds into the simulation directory ($TMPDIR) on the node.
Then a while loop copies all dirs called run_* back to the submission directory every half hour. The script named `kmc_master_parallel.py` is then executed with mpirun. `kmc_master_parallel.py` should look similar to the :ref:`quick-start example <quick>`:

.. code-block:: python

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


For further explanation, please see the `ASE documentation on parallel calculations <https://wiki.fysik.dtu.dk/ase/ase/parallel.html>`_.




