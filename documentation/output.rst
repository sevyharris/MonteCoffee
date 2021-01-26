.. _output:
.. index:: Code output


Code output
**************


Simulation log file
----------------------
:program:`MonteCoffee` generates a .txt log file with the naming convention kMClog_#date_#time_.txt. This is handled by the  module `base.Logging <api/NeighborKMC.base.html#module-NeighborKMC.base.logging>`_.

The first output contains all the dictionary of parameters, passed as the input :code:`parameters` when instantiating a :class:`NeighborKMC` object. This section may look as:

.. code-block:: python

    Mikkel Jorgensen (2015-2019)
    Noemi Bosio (since 2019)
    Elisabeth M. Dietze (since 2020)
    Chalmers University of Technology
    Goteborg, Sweden
    --------------------------------------------------------------------------------
    Simulation parameters
    pCO       :    2000.0
    pO2       :    1000.0
    T         :    800.0
    Name      :    COOx Simulation
    reverses  :    {0: 1, 2: 3, 4: 4, 5: 5}
    tend      :    1e-07
    Nsites    :    432
    --------------------------------------------------------------------------------
    No time acceleration used.

Next the log begins with printing the simulation step, the time the log-point was dumped, the simulation time, and the number of events fired, respectively:

.. code-block:: python

    Step          time[hr:min:s]             Sim time [s]          Events called
    1000             14:47:40             4.84677160781e-10        ['37', '0', '3', '0', '960', '1', '0']
    2000             14:47:41             7.10054219544e-10        ['50', '0', '4', '0', '1945', '2', '0']
    3000             14:47:42             8.66046768567e-10        ['59', '0', '4', '0', '2936', '2', '0']

The ordering of `Events Called` is determined by the order of which the parameter `events` are passed upon instantiating a :class:`NeighborKMC` object.

Results as .txt files
------------------------

:program:`MonteCoffee` generates its output as plain-text (.txt) file for basic output and as object in the HDF5 format for more advanced output, as determined by the method
`save_txt() in base.NeighborKMC <api/NeighborKMC.base.html#NeighborKMC.base.logging.NeighborKMCBase.save_txt>`_.

If desired, this method can be altered to output other information during runtime. The output files are generated using `numpy's savetxt() method <https://docs.scipy.org/doc/numpy/reference/generated/numpy.savetxt.html>`_, and can therefore be loaded using `np.loadtxt() <https://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html>`_.

The code outputs the following .txt files once prior to simulation:

    - **siteids.txt**: The index of each site, passed as the parameter :code:`ind` when instantiating a :class:`NeighborKMC.user_sites.Site` object. This is useful for storing an :class:`Ase.Atoms` object.
    
    - **stypes.txt**: The site types for each site, passed as the parameter :code:`stype` when instantiating a :class:`NeighborKMC.user_sites.Site` object.
    
During the simulation the following .txt files are updated with a frequency specified in :ref:`Options <options_sec>`:

    - **mcstep.txt**: The Monte Carlo step corresponding to each line in the other .txt files. This is directly dependent on the updating frequency.
    - **time.txt**: The simulation time in seconds for every logged Monte Carlo step.
    - **coverages.txt**: The coverages at each time-step for each lattice-site. The coverages are structured as :code:`coverages[mcstep][site_number]`.
    - **evs_exec.txt**: The total number of event-type executions. For example, to find the total number of executions of event number 0:
   
      .. code-block:: python
       
           evs_exec = np.loadtxt("evs_exec.txt")
           N1 = evs_exec[0]

  
    - **detail_site_event_evol.hdf5**: The detailed time evolution of the system, which saves every step. The file can be read in as follows:

      .. code-block:: python
      
           with h5py.File('detail_site_event_evol.hdf5','a') as f2:
               time_list = np.array(f2["time"]) 
               event_list = np.array(f2["event"]) 
               site_list = np.array(f2["site"]) 
               othersite_list = np.array(f2["othersite"]) 

     Working with this file, one has to be careful, because it can be rather big, despite the compression. In principle hdf5 can be read also with any other hdf5 supporting 
     programming language. 
 
For further information about analyzing output, see :ref:`analyzecoox` and :ref:`Calculating a turnover frequency <tof>`.
    
    

