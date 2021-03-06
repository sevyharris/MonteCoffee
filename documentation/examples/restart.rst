.. _tof:
.. index:: Restarting MonteCoffee

Restarting your simulation
*************************************

To collect sufficient statistic for the kinetic MonteCarlo simulation it may be
important to restart MonteCoffee after a certain time has elapsed, especially
on computer clusters with short queue times. 

Thus, the user part can be adjusted in the following way, first in the 
main-run file:

.. code-block:: python

    real_t_end = 10 #Real end time of simulation to restart in s
    # Instantiate simulator object, now including the simulation end time.
    sim = NeighborKMC(system=p, tend=tend,
                      real_t_end = real_t_end,
                      parameters=parameters,
                      events=events,
                      rev_events=reverse_events)

In the 'user_kmc.py' two new functions need to be defined, serialize and 
deserialize and the package pickle imported:

.. code-block:: python

   import pickle, os, time 

   def serialize(self,filename):
        """Ads the possibility to dump self object"""
        with open(filename, 'wb') as f:
            pickle.dump(self.__dict__,f)


    def deserialize(self,filename):
        """Reads the self object from the file"""
        with open(filename, 'rb') as f:
            self.__dict__ = pickle.load(f)

Additionally, the variable real_t_end has to be added to the __init__ of the simulation:

.. code-block:: python

      def __init__(self, system, tend, real_t_end = (96*60*60), parameters={}, events=[], rev_events={}):
        self.events = [ev(parameters) for ev in events]
        self.reverses = None # Set later
        self.real_t_end = real_t_end

The time module is used to follow the real time of the simulation. To use the real time
as second break condition of the simulation, it is included in the while-clause. At the 
end of the while-clause the self-object with the system state is dumped as pickle-file. 

.. code-block:: python

        start_time = time.time()  # save start time of simulation 
        if os.path.exists('data.pck'):   # if restart file exists, load self-object
            self.deserialize('data.pck')

        log.dump_point(self.stepNMC, self.t, self.evs_exec)
        while  time.time() < start_time + self.real_t_end and self.t < self.tend: 
            self.frm_step()

        self.serialize('data.pck') # dump self-object

Please notice: The time used here is the bare simulation time. Thus it must be reduced by
any pre-process time to initialize the system. 
