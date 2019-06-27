.. _options_sec:
.. index:: Options


Options
**************

The configurations and options in :program:`MonteCoffee` are specified in the file `kMC_options.cfg`, and are loaded by
`load_options() <api/NeighborKMC.base.html#NeighborKMC.base.kmc.NeighborKMCBase.load_options>`_.
The options contained in `kMC_options.cfg` are:

    - **nspecies** (int): how many species are included in the simulation. This integer decides how many coverages are calculated in `get_coverages() <api/NeighborKMC.base.html#NeighborKMC.base.kmc.NeighborKMCBase.get_coverages>`_.
    
    - **nninteractions** (int): describes the extent of the nearest neighbor interactions. Specifying 1 makes the local search for newly enabled events extent to second nearest neighbor. See `frm_update() <api/NeighborKMC.base.html#NeighborKMC.base.kmc.NeighborKMCBase.frm_update>`_.
    
    - **savesteps** (int): Number of Monte Carlo steps between saving the .txt files.
    
    - **logsteps** (int): Number of steps to skip when saving .txt files and printing the log.
    
    - **tinfinity** (float): The time to consider infinite. This time is set for the occourence time of impossible events.
    
    - **savecovs** (bool): If True, the site-occupations (coverages) are saved together with the other .txt files.
    
    - **verbose** (bool): If True, prints verbose information.
    
Options solely related relates to the acceleration of kMC simulations (See :ref:`accelerating <accelerating>`):
    
    - **delta** (float): The reversibility tolerance for quasi-equilibrated events. This determines the criterion for whether reactions are quasi-equilibrated.
    
    - **nf** (int): The average number of steps to take for each quasi-equilibrated event in the superbasin.
    
    - **ns** (int): The frequency of adjusting barriers for quasi-equilibrated events.
    
    - **ne** (int): How many events back to remember when adjusting barriers of quasi-equilibrated events.
    
    - **usekavg** (bool): If True, the average rate-constant of non-equilibrated events in the superbasin is used to scale the barriers of the quasi-equilibrated events, and not the actual rates. Setting `usekavg` to True provides a more stable time-step, but can be overly conservative.
