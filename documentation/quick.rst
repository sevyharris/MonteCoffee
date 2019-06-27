.. _quick:
.. index:: Quick start


Quick start
**********************************
This tutorial is supposed to provide a minimal working example.
For a full in-depth tutorial demonstrating CO oxidation over a Pt nanoparticle,
see :ref:`The CO oxidation tutorial <coox>`. To keep the clarity high, all energy barriers are assumed constant,
only two types of reaction events are implemented, and neighbors are based on distances in an :class:`ASE.Atoms` object.
Moreover, everything is done in one python script.

**Step 1. Implement the two event-types**

In this step we import the template class :class:`NeighborKMC.base.events.EventBase`, and derive two classes from this to
enable two different types of reactions. All defined events must implement three methods: 

    - `possible(self, system, site, other_site) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.possible>`_
    - `get_rate(self, system, site, other_site) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.get_rate>`_
    - `do_event(self, system, site, other_site) <api/NeighborKMC.base.html#NeighborKMC.base.events.EventBase.do_event>`_

The first event class is an adsorption, that is possible if the pair of neighbor sites are empty, has a rate linearly dependent on the pressure, and 
if executed, it covers the sites with species 1:

>>> from base.events import EventBase
>>>
>>> class Adsorption(EventBase):
>>>
>>>     def __init__(self, params):
>>>         EventBase.__init__(self, params)
>>>
>>>     def possible(self, system, site, other_site):
>>>
>>>         if system.sites[site].covered == 0 and system.sites[other_site].covered == 0:
>>>             return True
>>>         else:
>>>             return False
>>>
>>>     def get_rate(self, system, site, other_site):
>>>         R = 1000. * self.params["pA"]
>>>         return self.alpha * R  # alpha important for temporal acceleration.
>>>
>>>     def do_event(self, system, site, other_site):
>>>         # Cover the two sites with species 1
>>>         system.sites[site].covered = 1
>>>         system.sites[other_site].covered = 1


Now we define the reverse desorption event with a constant rate:

>>> class Desorption(EventBase):
>>>
>>>     def __init__(self, params):
>>>         EventBase.__init__(self, params)
>>>
>>>     def possible(self, system, site, other_site):
>>>
>>>         if system.sites[site].covered == 1 and system.sites[other_site].covered == 1:
>>>             return True
>>>         else:
>>>             return False
>>>
>>>     def get_rate(self, system, site, other_site):
>>>         R = 101.
>>>         return self.alpha * R  # alpha important for temporal acceleration.
>>>
>>>     def do_event(self, system, site, other_site):
>>>         # empty the sites:
>>>         system.sites[site].covered = 0
>>>         system.sites[other_site].covered = 0

Now we will store **references** to the classes in a list, and make a `dict` that shows that these two events are reverses:

>>> events = [Adsorption, Desorption]
>>> # Specify what events are eac others' reverse.
>>> reverse_events = {0: 1}

Specifying which events are reverses help accelerating kMC simulations, as described in the tutorial on :ref:`Accelerating kMC <accelerating>`.

**Step 2. Define sites**

In this step, the sites are defined from an :class:`ASE.Atoms` object. Here we create one site for each atom in
a 10x10 fcc(111) surface, all with the same site-type :code:`stype`:

>>> from ase.build import fcc111
>>> from user_sites import Site
>>> 
>>> a0 = 4.00  # Lattice Parameter (not related to DFT!)
>>> atoms = fcc111("Pt", size=(10, 10, 1), a=a0)
>>> sites = []
>>> # Define a site for each atom that is empty with no pre-defined neighbors:
>>> for i in range(len(atoms)):
>>>     sites.append(Site(stype=0, covered=0, ind=[i]))

Now we have a list of empty sites, which are used to instantiate a system.

**Step 3. Instantiate system and neighborlists**

Here, the system is created and the sites are connected by calcualting a neighborlist. In this example,
the `set_neighbors() <api/NeighborKMC.html#NeighborKMC.user_system.System.set_neighbors>`_ method is used, which assigns sites that are separated by no more than one nearest neighbor distance:

>>> import numpy as np
>>> from user_system import System
>>> p = System(atoms=atoms, sites=sites)
>>> Ncutoff = a0 / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
>>> p.set_neighbors(Ncutoff, pbc=True)

:code:`pbc=True` turns on :ref:`periodic boundary conditions <pbc>`.

**Step 4. Instantiate a NeighborKMC object and run**

Now we are ready to instantiate a :class:`NeighborKMC.user_kmc.NeighborKMC` object, which is connecting the ingredients created in the previous step.  But first we create a `dict` containing all the parameters passed onto the events to calculate rates:

>>> from user_kmc import NeighborKMC
>>> parameters = {"pA": 100., "Name": "Quickstart simulation", "reverses ": reverse_events}
>>> sim = NeighborKMC(system=p, 
>>>                   tend=1.0, # end after 1.0 s.
>>>                   parameters=parameters, # parameters for event rate-constants.
>>>                   events=events, # the list of events
>>>                   rev_events=reverse_events) # the dict of reverse events
 
 
Now we can run the simulation by invoking

>>> sim.run_kmc()

Then it is just to have a cup of coffee and wait.

**Afterthoughts**

While this example shows how simple it can be to run a simulation, in more complex examples it is useful to store the information separate files:

    - `user_events.py <api/NeighborKMC.html#module-NeighborKMC.user_events>`_ can be used to store the event-types.
    - `user_energy.py <api/NeighborKMC.html#module-NeighborKMC.user_energy>`_ can be used to store functions for obtaining energies used to calculate event rate constants.
    - `user_entropy.py <api/NeighborKMC.html#module-NeighborKMC.user_entropy>`_ can be used to store entropy calculation functions.
    - `user_constants.py <api/NeighborKMC.html#module-NeighborKMC.user_constants>`_ can be used to store global and physical constants.
    
    
A full example following these guidelines is shown in the tutorial :ref:`CO oxidation on a Pt nanoparticle <coox>`.
