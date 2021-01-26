.. _parallel:
.. index:: Parallelisation

Parallel simulations
*************************************
Kinetic MonteCarlo simulations are stochastic in nature, making it reasonable to perform multiple identically prepared simulations
to assess the convergence. In the following, we assume that `ASE <https://wiki.fysik.dtu.dk/ase/>`_ and `MPI4PY <https://pypi.org/project/mpi4py/>`_ are installed.

To submit a simulation, defined in a file `kmc_master_parallel.py`, in an environment that implements `SLURM <https://slurm.schedmd.com/>`_, one can submit a bash script as

.. code-block:: bash

    #!/usr/bin/env bash
    #SBATCH -N 1
    #SBATCH -n 10
    #SBATCH -p PROJNAME
    #SBATCH -A ACCNO
    #SBATCH -o out.txt
    #SBATCH -e err.txt
    #SBATCH -t 96:15:00
    #SBATCH -J project_dir/submit_dir
    #SBATCH --mail-user=USER@UNI.se
    #SBATCH --mail-type=END

    cp *.py $TMPDIR
    cp kMC_options.cfg $TMPDIR
    cp -r base $TMPDIR

    cd $TMPDIR

    while sleep 1800; do
        # This will be executed once per every 3600 seconds
        rsync -a $TMPDIR/* $SLURM_SUBMIT_DIR
    done &     # The &-sign after the done-keyword places 
               # the while-loop in a sub-shell in the background
    LOOPPID=$! # Save the PID of the subshell running the loop

    mpirun -np 10 python kmc_master_parallel.py 

    # Copy the files back after run:
    rsync -a $TMPDIR/* $SLURM_SUBMIT_DIR 
    kill $LOOPPID

    
    
This script copies all python files in MonteCoffee to a compute node, and cds into the simulation directory ($TMPDIR) on the node.
Then a while loop copies all dirs called run_* back to the submission directory every half hour. The script named `kmc_master_parallel.py` is then executed with mpirun. The script `kmc_master_parallel.py` can at first be very similar to the :ref:`quick-start example <quick>`. For simplification the detailed definitions as in :ref:`quick-start <quick>` are not given her, so that one can see the differences. Mainly they are in the first part of the code. 

.. code-block:: python

    import os
    from ase.build import fcc111
    from ase.parallel import MPI4PY
    from base.sites import SiteBase
    import numpy as np
    from base.system import SystemBase
    from base.logging import Log
    from base.kmc import NeighborKMCBase

    world = MPI4PY()
    rank = world.rank # What number simulation copy am I?
    size = world.size # How many total simulation copies?
    rundir = "run_"+str(rank)  # Name of the dir I create
    
    os.mkdir(rundir) # Create dir
    os.system("cp kMC_options.cfg "+rundir) # copy kMC_options.cfg here
    os.chdir(rundir)

    a0 = 4.00  # Lattice Parameter (not related to DFT!)
    Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff

    atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
    sites = []

    for i, indic in enumerate(atoms):
        sites.append(Site(stype=0, covered=0, ind=[i]))

    events = [Adsorption, Desorption]
    
    p = SystemBase(atoms=atoms, sites=sites)
    p.set_neighbors(Ncutoff)

    parameters = {"pA": 10., "Name": "Parallel Simulation"}

    sim = simple_NKMC(system=p,
                      tend=10.,
                      parameters=parameters, 
                      events=events)
                      
    sim.run_kmc()


For further explanation about using MPI4PY within ASE, please see the `ASE documentation on parallel calculations <https://wiki.fysik.dtu.dk/ase/ase/parallel.html>`_. 

In general, it can be useful to assign a large :code:`tend` and let the bash-script runtime determine the end of simulation. Because the code itself writes out log-files regularly, one will not loose any informations by letting the script runtime determine the end of the simulation.


