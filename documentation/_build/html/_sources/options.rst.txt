.. _options_sec:
.. index:: Options


Options
**************

The configurations and options in :program:`MonteCoffee` are specified in the file `kMC_options.cfg`, and are loaded by
`load_options() <api/NeighborKMC.base.html#NeighborKMC.base.kmc.NeighborKMCBase.load_options>`_.
The options contained in `kMC_options.cfg` are:

    - **nspecies** (int): how many species are included in the simulation. This integer decides how many coverages are written to the log.

    - **nninteractions** (int): describes the extent of the nearest neighbor interactions. Specifying 1 makes the local search for newly enabled events extent to second nearest neighbor.

    - **savesteps** (int): Number of Monte Carlo steps between saving the .txt files.
    
    - **logsteps** (int): Number of steps to skip when saving .txt files and printing the log. Setting it to 1000 makes the log output and code save every 1000th simulation step.
    
    - **tinfinity** (float): The time to consider infinite. This time is set for the occurrence time of impossible events.
    
    - **savecovs** (bool): If True, the site-occupations (coverages) are saved periodically to coverages.txt.
    
    - **verbose** (bool): If True, prints verbose information.

    - **write_atoms** (bool): If True, write out traj files of coverages (may have to be customized depending on the system).
    
Options solely related to the acceleration of kMC simulations (See :ref:`accelerating <accelerating>`):

    - **use_scaling_algorithm** (str): If None or False, no time acceleration will be used. In the current version the following three schemes are available:

        - :code:`scale_constant` -- scaling with a constant factor
        - :code:`scale_rate` -- scaling based on the reaction rates
        - :code:`scale_rate_constant` -- scaling based on the reaction rate constants
    
    - **delta** (float): The reversibility tolerance for quasi-equilibrated events. This determines the criterion for whether reactions are quasi-equilibrated. Setting delta=0.25 implies that there must be no more than 25 percent difference in rates of a forward and backward reaction to deem it quasi-equilibrated.
    
    - **nf** (int): A factor that separates fast and slow events. nf = 1000 means that quasi-equilibrated events are slowed down to maximally 1000 times the non-equilibrated events.
    
    - **ns** (int): The frequency of adjusting barriers for quasi-equilibrated events.
    
    - **ne** (int): The number of events executed for a reaction in forward and backward direction for the reaction to be considered equilibrated. 
