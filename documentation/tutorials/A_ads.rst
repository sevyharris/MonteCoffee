.. _A_ads:
.. index:: A adsorption

Adsorption of atomic species A on surface
******************************************

In this first tutorial the first steps to run and set-up a :program:`MonteCoffee` simulation are shown, which require only a view steps. 
.. The entire tutorial is shown in `test.py <../api/NeighborKMC.html#module-NeighborKMC.test>`_ and the references to the other modules mentioned herein.

In this guide the simple adsorption/desorption of an atomic species on a plain surface is demonstrated and the results compaired to the solution of a mean-field model. The reaction to simulate is:

.. math::

   A + * \longleftrightarrow A^*

and the time evolution of the coverage of species A according to the mean-field model (see e.g.: `Fichthorn and Weinberg <http://aip.scitation.org/doi/10.1063/1.461138>`_)

.. math::

   \theta(t) = \frac{r_A}{r_A+r_D}(1-e^{-(r_A+r_D)*t})

with :math:`\theta` the coverage of species A, :math:`t` the time and :math:`r_{A,D}` the rate of adsorption and desorption respectively. 

Before the simulation, constants, reaction sites and system as well as the events have to defined which will be shown in the next steps.  

First, various  parameters must be set and which can be stort in a parameter-dictionairy:

.. code-block:: python

     tend = 10.  # End time of simulation (s)
     latt_param = 4.00  # Lattice Parameter (not related to DFT!)
     Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff


Define sites and system
----------------------------
One site is defined for each surface atom using an :class:`Ase.Atoms` object.
We start with an empty list of sites and a fcc(100) surface in `ASE`:

.. code-block:: python

     from ase.build import fcc(100)
     from user_sites import Site

     surface = fcc100("Pt", a=latt_param, size=(10,10,1))
     sites = []

Now we can create a site, free of adsorbates, for each surface atom with a :code:`stype` . Here all surface atoms are of the same kind.:

.. code-block:: python

     # Create a site for each surface-atom:
     for i in range(len(atoms)):
         sites.append(Site(stype=0,
                           covered=0, ind=[i]))

Here, the block :code:`ind=[i]` stores the index of the atom in the :class:`ASE.Atoms` object on the :class:`NeighborKMC.user_sites.Site` object.

Finally, we need to define neighborlists. It is simplest to define this according to the nearest neighbor distances:

.. code-block:: python

     # Set the neighbor list for each site using distances.

     for i, s in enumerate(sites):
         for j, sother in enumerate(sites):
             # Length of distance vector:
             dcur = self.atoms.getdistance(s.ind[0], sother.ind[0], mic = pbc)

             # If the site is a neighbor:
             if dcur < Ncutoff and j != i:
                 s.neighbors.append(j)

 
Now the :class:`NeighborKMC.user_system.System` object can be defined from the collection of sites:

.. code-block:: python

     from user_system import System
     p = System(atoms=atoms, # store ASE.Atoms as well
                sites=sites)

Define events
--------------
Here event-types are defined, which are stored in `user_events.py <../api/NeighborKMC.html#module-NeighborKMC.user_events>`_.
For each possible type of event, a class is derived from :class:`NeighborKMC.base.events.EventBase`. In this case, we need to define two different events, the adsorption of species A, and correspondingly the desorption. 

First we import the necessary functions, classes, and constants:

.. code-block:: python

     from base.events import EventBase

Now we derive a class to contain the event:

.. code-block:: python

     class AAdsEvent(EventBase):
         def __init__(self, params):
             EventBase.__init__(self, params)

The constructor :code:`__init__(self,params)` is attaches relevant parameters to the object. 
We need a function `possible(self,system, site, other_site)` that returns True if the event is possible on the current site-pair. For single atom adsorption it does not matter if the other_site is covered or not. Thus we are only interested in the site itself.

.. code-block:: python

         def possible(self, system, site, other_site):
             # If site is uncovered 
             if (system.sites[site].covered == 0):
                 return True
             else:
                 return False

Thus, for the event to be possible, the site needs to be empty.
Now we also need to define a function :code:`get_rate(self, system, i_site, other_site)` that returns the rate constant. To keep this as simple as possible, the rate constant is chosen to be :math:`R=1`.

.. code-block:: python

        def get_rate(self, system, i_site, other_site):
            R = 1.
            return R


Finally each event requires a method :code:`do_event(self,system, site, other_site)` to perform modifications to the site-occupations when fired:

.. code-block:: python

        def do_event(self, system, site, other_site):
            system.sites[site].covered = 1

In this case, up on adsorption the site is covered with the species A, represented by the number 1 within the code. 

To take care of the correct time evolution of the :program:`MonteCoffee` we introduce an aditional block with returns if either neighboring sites are involved or not. Here no neighbouring sites are involved, thus we :code:`return False`. 

.. code-block:: python 

        def get_involve_other(self):
            return False 

To know which events are each others reverse we store:

.. code-block:: python

     reverse_events = {0: 1}
     events = [AAdsEvent, ADesEvent]

Here event 0 has a reverse event 1.
The numbering of events is determined by the order in the list :code:`events` defined here.

Define and run simulation
-----------------------------

Now the simulation object :class:`NeighborKMC.user_kmc.NeighborKMC` can be defined and the simulation performed:

.. code-block:: python

   parameters = { "Name": "A ads/des Simulation",
                  "reverses ": reverse_events}

     # Instantiate simulator object.
     sim = NeighborKMC(system=p, tend=tend,
                       parameters=parameters,
                       events=events,
                       rev_events=reverse_events)
     result = sim.run_kmc()
     print("Simulation end time reached ! ! !")




.. _analyzecoox:

Analyze results
----------------------------
The results are analyzed by reading in the :ref:`code output <output>`. Here, we would like to calculate the A coverage as a function of time for the entire system:

.. code-block:: python

     import numpy as np
     time = np.loadtxt("time.txt")
     covs = np.loadtxt("coverages.txt")
     Nsites = float(len(covs[0]))
     cov_A = [sum([1 for val in covs[i] if val == 1]) / Nsites for i in range(len(covs))]

This can be plotted as done in the following example with :code:`matplotlib`

.. code-block:: python

     import matplotlib.pyplot as plt
     plt.plot(time, cov_A, '-k')
     plt.xlabel("Time [s]")
     plt.ylabel("Coverage")
     plt.savefig('coverage_spec_A.pdf')

To compare the effect of the used simulation surface on the result and also compare to the result of the mean-field model in the following a plot is shown with surface sizes of (5x5), (10x10) and (100x10) corresponding to 25, 100 and 1000 surface sites respectively. Also shown is the result of a very simple first-reaction kinetic MonteCarlo code with 1000 sites. 

.. image:: ../images/compare_MF_kMC.pdf 

If an increase in the number of sites is not possible, it is recommanded that multiple identically prepared simulations are performed.
(see  tutorials :ref:`Parallel simulations <parallel>` and :ref:`calculating turnover frequencies <tof>`).





