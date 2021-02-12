.. _quick:
.. index:: Quick start


Quick start
**********************************
This tutorial provides a minimal working example of how to run :program:`MonteCoffee`. For simplicity this example puts all definitions in one file. This structure is different from the following tutorials.  

For this very first, simple example, all free energy barriers are assumed constant. Two different events are implemented, and neighbors are calculated from inter-atomic distances in an :class:`ASE.Atoms` object. For simplicity, this guide shows how to setup a simulation in a single python script. In principle it is similar to the second tutorial on diatomic adsorption and desorption (see :ref:`b2ads_tut`).

**Step 1. Implement the two event-types**

In this step we import the template class :class:`NeighborKMC.base.events.EventBase`, and derive two classes from this to
enable two different types of reactions. All defined events must implement four methods: 

    - `possible(self, system, site, other_site) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.possible>`_: returns True if an event is possible.
    - `get_rate(self, system, site, other_site) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.get_rate>`_: returns the rate constant of the event.
    - `do_event(self, system, site, other_site) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.do_event>`_: executes the event by modifying site-occupations.
    - `get_involve_other(self) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.get_involve_other>`_: defines if neighboring sites are important for event.


The first event class is an adsorption, that is possible if the pair of neighbor sites is empty.
Each eventclass can be given a name which is only used in the output.
The adsorption class has a rate-constant that is linearly dependent on the pressure, and
if executed, it covers the sites with species 1:

.. code-block:: python

    from base.events import EventBase

    class Adsorption(EventBase):

        def __init__(self, params):
            EventBase.__init__(self, params, name='Adsorption')

        def possible(self, system, site, other_site):

            if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
                return True
            else:
                return False

        def get_rate(self, system, site, other_site):
            R = 1. * self.params["pA"]
            return R  

        def do_event(self, system, site, other_site):
            # Cover the two sites with species 1
            system.sites[site].covered = 1
            system.sites[other_site].covered = 1

        def get_involve_other(self):
            return True


Now we define the reverse desorption-event with a constant rate:

.. code-block:: python

    class Desorption(EventBase):

        def __init__(self, params):
            EventBase.__init__(self, params, name='Desorption')

        def possible(self, system, site, other_site):

            if system.sites[site].covered == 1 and system.sites[other_site].covered == 1:
                return True
            else:
                return False

        def get_rate(self, system, site, other_site):
            R = 1.
            return R  

        def do_event(self, system, site, other_site):
            # empty the sites:
            system.sites[site].covered = 0
            system.sites[other_site].covered = 0

        def get_involve_other(self):
            return True

Now we will store **references** to the classes in a list: 

.. code-block:: python

    events = [Adsorption, Desorption]

How to accelerate the kMC simulations is described in the tutorial on :ref:`Accelerating kMC <accelerating>` and not used in this simplest case.

**Step 2. Define sites**

In this step, the sites are defined from an :class:`ASE.Atoms` object. We create one site for each atom in
a 10x10 fcc(111) surface, all with the same site-type :code:`stype=0` and without any covering species :code:`covered=0`:

.. code-block:: python

    from ase.build import fcc111
    from base.sites import SiteBase

    a0 = 4.00  # Lattice Parameter (not related to DFT!)
    atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
    sites = []
    # Define a site for each atom that is empty with no pre-defined neighbors:
    for i in range(len(atoms)):
        sites.append(SiteBase(stype=0, covered=0, ind=i))

Now we have a list of empty sites, which are used to instantiate a system.

**Step 3. Instantiate system and neighborlists**

Here, the system is created and the sites are connected by calculating a neighborlist. In this example we assign neighbors within one nearest-neighbor distance:

.. code-block:: python

    import numpy as np
    from base.system import SystemBase

    p = SystemBase(atoms=atoms, sites=sites)
    Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff

    for i, s in enumerate(sites):
        for j, sother in enumerate(sites):
            dcur = atoms.get_distance(s.ind, sother.ind, mic=True)
            if dcur < Ncutoff and j != i:
                s.neighbors.append(j)

:code:`mic=True` uses periodic boundary conditions to imply an infinite surface. 

**Step 4. Instantiate a NeighborKMC object and run**

Now we are ready to instantiate a :class:`NeighborKMC.NeighborKMCBase` object, which is connecting the ingredients created in the previous step. The main part of the kinetic Monte Carlo procedure is in the :class:`NeighborKMC.NeighborKMCBase`, but some details and logging should be defined by the user. That is done here in the class: :class:`simple_NKMC`.

.. code-block:: python

    from base.logging import Log
    from base.kmc import NeighborKMCBase

    class simple_NKMC(NeighborKMCBase):

        # First we initialize the kMC simulation and load the parameters
        def __init__(self, system, tend, parameters={}, events=[]):
            self.events = [ev(parameters) for ev in events]
            NeighborKMCBase.__init__(self, system=system,
                                     tend=tend, parameters=parameters)
    
       # We also define a run_kmc, which runs the actual kMC simulation. Also the logging is defined here. 
        def run_kmc(self):     
            logparams = {}
            logparams.update(self.parameters)
            logparams.update({"tend": self.tend,
                               "Nsites": self.system.Nsites,
                               "Number of events": len(self.events),
                               "Number of site-types (stypes)": len(list(set([m.stype for m in self. system.sites])))
                               })
            log = Log(logparams)
    
            stepN_CNT = 0  # Parameter to count LogSteps threshold
            stepNMC = 0    # Parameter to count the number of executed kMC steps
    
            while self.t < self.tend:  
                self.frm_step()         # Execute a kMC step
    
                if stepN_CNT >= self.LogSteps:       # Only for Logging purposes 
                    print("Time : ", self.t, "\t Covs :", self.system.get_coverages(self.Nspecies))
                    log.dump_point(stepNMC, self.t, self.evs_exec)
                    stepN_CNT = 0
    
                stepN_CNT += 1
                stepNMC += 1
 

So now after defining the :class:`simple_NKMC` object, the parameters can be loaded:

.. code-block:: python

    parameters = {"pA": 10., "Name": "Quickstart simulation"}
    sim = simple_NKMC(system=p,
                      tend=10.0, # end after 10.s.
                      parameters=parameters, # parameters for event rate-constants.
                      events=events) # the list of events
 
And finally, we can now run the simulation by invoking:

.. code-block:: python

    sim.run_kmc()

Then it is just to have a cup of coffee and wait.

**Afterthoughts**

While this example shows how simple it can be to run a simulation, in all following tutorials and examples the details of the simulations are stored in so-called
user-files:

    - `user_kmc.py <api/NeighborKMC.html#module-NeighborKMC.user_kmc>`_ can be used to customize the kMC routine, especially run_kmc().
    - `user_events.py <api/NeighborKMC.html#module-NeighborKMC.user_events>`_ can be used to store the event-types.
    - `user_energy.py <api/NeighborKMC.html#module-NeighborKMC.user_energy>`_ can be used to store functions for obtaining energies used to calculate event rate constants.
    - `user_entropy.py <api/NeighborKMC.html#module-NeighborKMC.user_entropy>`_ can be used to store entropy calculation functions.
    - `user_constants.py <api/NeighborKMC.html#module-NeighborKMC.user_constants>`_ can be used to store global and physical constants.
    
