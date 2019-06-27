.. _coox:
.. index:: CO oxidation

CO oxidation on a Pt nanoparticle
*************************************

To begin running :program:`MonteCoffee` simulations, there are only a few required steps.
The entire procedure is shown in `test.py <../api/NeighborKMC.html#module-NeighborKMC.test>`_ with references to the other modules mentioned herein.

The following guide shows how to perform CO oxidation simulations on a Pt truncated octahedron with an energy landscape based on conventional coordination numbers.
By going through this, all information that is relevant to performing simulations on more complex reaction networks should be revealed.

In this guide, species 0 is free sites, 1 is CO and 2 refers to O.

Define constants
----------------------
First any global constants can be defined in `user_constants.py <../_modules/NeighborKMC/user_constants.html>`_.
Next parameters are set and stored in a dictionary as:

>>> T = 800.  # Temperature
>>> pCO = 2E3  # CO pressure
>>> pO2 = 1E3  # O2 pressure
>>> tend = 1E-3  # End time of simulation (s)
>>> parameters = {"pCO": pCO, "pO2": pO2, "T": T,
>>>               "Name": "COOx Simulation",
>>>               "reverses ": reverse_events}
>>> a = 4.00  # Lattice Parameter (not related to DFT!)
>>> Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff


Define sites and system
----------------------------
In this step, the sites are defined. For example, one site can be defined for each surface atom using an :class:`Ase.Atoms` object:
First, we start with an empty list of sites and a nanoparticle in `ASE`:

>>> from ase.cluster import Octahedron
>>> from user_sites import Site
>>>
>>> sites = []
>>> particle = Octahedron("Pt", length=12, cutoff=3, latticeconstant=a)

Now we define one site for each surface atom in particle by first finding the coordination number:

>>> # First find CN of each atom:
>>> CNS = np.zeros(len(atoms))
>>> for i, at in enumerate(atoms):
>>>     pcur = at.position
>>>     dp = np.sqrt([(p[0] - pcur[0]) ** 2. + (p[1] - pcur[1]) ** 2. +
>>>                   (p[2] - pcur[2]) ** 2. for p in atoms.positions])
>>>     CNS[i] = len([val for val in dp if 0. < val < a / np.sqrt(2) + 0.01])
>>>
>>> # Define surface atoms as non-bulk.
>>> surface_atom_ids = [i for i in range(len(CNS)) if CNS[i] < 12]

In this example, let each coordination number corresponds to a unique site-type (`stype`):

>>> stypes = {}  # site-types: (111),(100),edge,corner.
>>> for i, k in enumerate(sorted(list(set(CNS)))):
>>>     stypes[k] = i

Now we can create a site, free of adsorbates, for each surface atom with a site-type that corresponds to the coordination number:

>>> # Create a site for each surface-atom:
>>> for i, indic in enumerate(surface_atom_ids):
>>>     sites.append(Site(stype=stypes[CNS[indic]],
>>>                       covered=0, ind=[indic]))

Here, the block :code:`ind=[indic]` stores the index of the atom in the :class:`ASE.Atoms` object on the :class:`NeighborKMC.user_sites.Site` object.


Finally, we need to define neighborlists. It is simplest to define this according to the nearest neighbor distances:

>>> # Set the neighbor list for each site using distances.
>>>
>>> for i, s in enumerate(sites):
>>>     pcur = atoms[s.ind[0]].position # position of site
>>>     for j, sother in enumerate(sites):
>>>         pother = atoms[sother.ind[0]].position # position of potential neighbor site:
>>>         # Length of distance vector:
>>>         dpabs = np.sqrt((pother[0] - pcur[0]) ** 2. +
>>>                         (pother[1] - pcur[1]) ** 2. +
>>>                         (pother[2] - pcur[2]) ** 2.)
>>>
>>>         # If the site is a neighbor:
>>>         if dpabs < Ncutoff and j != i:
>>>             s.neighbors.append(j)

 
Now the :class:`NeighborKMC.user_system.System` object can be defined from the collection of sites:

>>> from user_system import System
>>> p = System(atoms=atoms, # store ASE.Atoms as well 
>>>            sites=sites)


Define reaction energies and entropies
--------------------------------------------
In this step, the reaction energies, or methods to calculate these, are defined in `user_energy.py <../api/NeighborKMC.html#module-NeighborKMC.user_energy>`_.
**In principle, one may skip this section** and simply define reaction energy barriers directly in :ref:`define events <defeventsquick>`, however, we believe this
step is good for keeping an overview of the coding of the energy landscape.

In this example from  `user_energy.py <../api/NeighborKMC.html#module-NeighborKMC.user_energy>`_, the adsorption energies of CO and O are stored as lists (functions of coordination number), the reaction energy barrier as a function `get_Ea(ECO, EO)`, and diffusion barriers as constants:

>>> EadsCO = [1.36 + 0.25 * (9 - CN) for CN in [6, 7, 8, 9]]
>>> EadsO = [0.97 + 0.2 * (9 - CN) for CN in [6, 7, 8, 9]]
>>>
>>> EdiffCO = 0.046 # CO diffusion barrier 
>>> EdiffO = 0.5 # O diffusion barrier
>>>
>>> def get_Ea(ECO, EO):
>>>    dEO = EO - EadsO[-1]  # Oxygen energy relative to uncovered Pt(111)
>>>    dECO = ECO - EadsCO[-1]  # CO energy relative to uncovered Pt(111)
>>>    dETS = 0.824 * (dEO + dECO)  # How much larger is the energy of CO and O wrt Pt(111)
>>>    Ea = 1.08 + dETS - dECO - dEO  # Translate the barriers relative to Pt(111)
>>>    return Ea

Repulsive adsorbate-adsorbate interactions are also defined as a method in  `user_energy.py <../api/NeighborKMC.html#module-NeighborKMC.user_energy>`_:

>>> def get_repulsion(cov_self, cov_NN, stype):
>>>     repulsion = 0.
>>>     ECOCO = 0.19  # 0.38 # How CO affects CO
>>>     EOO = 0.32  # How O affects O - double since it is called from get barrier of O2
>>>     ECOO = 0.3  # How CO affects O
>>>     EOCO = 0.3  # How O affects CO
>>>     HInttwo = [[0., 0., 0.], [0., ECOCO, EOCO],
>>>               [0., ECOO, EOO]]  # Two body interaction Hamiltonian 3x3 beacuse 0 = empty.
>>>     for j in cov_NN:  # For each covered Neighbor, give a repulsion:
>>>         repulsion += HInttwo[cov_self][j]
>>>
>>>     return repulsion

Now entropies are stored in `user_entropy.py <../api/NeighborKMC.html#module-NeighborKMC.user_entropy>`_, where the entropy is defined for gas-phase CO and oxygen, as well as a method to calculate harmonic adsorbate entropy. For brevity, please see the module `user_entropy.py <api/NeighborKMC.html#module-NeighborKMC.user_entropy>`_ for definition of the entropy functions.

.. _defeventsquick:

Define events
--------------
Here event-types are defined, which are stored in `user_events.py <../api/NeighborKMC.html#module-NeighborKMC.user_events>`_.
For each possible type of event, a class is derived from :class:`NeighborKMC.base.events.EventBase`, which requires defining 3 different methods.
Take the example of an event where CO+O forms CO2. This event is defined in `user_events.py <../api/NeighborKMC.html#module-NeighborKMC.user_events>`_ as follows.

First we import the necessary functions, classes, and constants:

>>> from base.events import EventBase
>>> from user_entropy import get_entropy_CO, get_entropy_O2, get_entropy_ads, get_Zvib
>>> from user_constants import mCO, mO2, Asite, modes_COads, modes_Oads, kB, eV2J, s0CO, s0O, h
>>> from user_energy import EadsCO, EadsO, get_Ea, get_repulsion, EdiffCO, EdiffO

Now we derive a class to contain the event:

>>> class COOxEvent(EventBase):
>>>     def __init__(self, params):
>>>         self.Zratio = (get_Zvib(params["T"], modes_COads) *
>>>                       get_Zvib(params["T"], modes_Oads)) ** 0.66
>>>         EventBase.__init__(self, params)

The constructor :code:`__init__(self,params)` is there to attach constants to the object, and :code:`self.Zratio` is the ratio
between the partition functions in the initial state and transition state, used to calculate the rate. Now we need to write a function
that returns True if the event is possible on the current site:

>>>     def possible(self, system, site, other_site):
>>>         # If site is covered with CO and other site free
>>>         if (system.sites[site].covered == 1 and
>>>                system.sites[other_site].covered == 2):
>>>             return True
>>>         else:
>>>             return False

Here we say that we need the site to be covered by 1 (CO) and the neighbor site by 2 (O) for the event to be possible.
Now we also need to define a function :code:`get_rate(self, system, i_site, other_site)` that returns the rate constant:

>>>    def get_rate(self, system, i_site, other_site):
>>>        # Find the adsorption energy for the site-type
>>>        stype = system.sites[i_site].stype
>>>        stype_other = system.sites[other_site].stype
>>>        ECO = EadsCO[stype]
>>>        EO = EadsO[stype_other]
>>>        # Find the Nearest neighbor repulsion
>>>        Ncovs = [system.sites[n].covered for n in
>>>                 system.neighbors[i_site]]
>>>        Nothercovs = [system.sites[n].covered for n
>>>                      in system.neighbors[other_site]]
>>>        ECO -= get_repulsion(1, Ncovs, stype)
>>>        EO -= get_repulsion(2, Nothercovs, stype_other)
>>>        Ea = max(0., get_Ea(ECO, EO)) # No negative energy barriers
>>>
>>>        return self.alpha * self.Zratio * np.exp(-Ea /
>>>                                                 (kB * self.params['T'])) * kB * self.params['T'] / h

Here, the site-types are used to get the adsorption energies, and the repulsions are added to the adsorption energies. Then a call is made to :code:`get_Ea(ECO, EO)` to obtain the reaciton energy barrier, and the rate is multiplied by the ratio of partition functions :code:`self.Zratio` to account for entropy losses. **It is important to multiply rate constants with** :code:`self.alpha` **if this event is supposed to be** :ref:`accelerated <accelerating>`. This is because :code:`self.alpha` is the slowing-down factor that is adjusted during simulation.

Finally each event requires a method :code:`do_event(self,system, site, other_site)` to perform modifications to the site-occupations when fired:

>>>    def do_event(self, system, site, other_site):
>>>        system.sites[site].covered = 0
>>>        system.sites[other_site].covered = 0

In this case, the two sites containing CO and O are simply emptied. Now, assume we have defined an event for each type of reaction desired:

    - (0) :class:`NeighborKMC.user_events.COAdsEvent` for CO adsorption.
    - (1) :class:`NeighborKMC.user_events.CODesEvent` for CO desorption.
    - (2) :class:`NeighborKMC.user_events.OAdsEvent` for O2 dissociative adsorption.
    - (3) :class:`NeighborKMC.user_events.ODesEvent` for O2 desorption.
    - (4) :class:`NeighborKMC.user_events.CODiffEvent` for CO diffusion.
    - (5) :class:`NeighborKMC.user_events.ODiffEvent` for O diffusion.
    - (6) :class:`NeighborKMC.user_events.COOxEvent` for CO+O -> CO2.

To accelerate the simulation we need to specify which events are each others inverse and store the event-class references in a list:

>>> reverse_events = {0: 1, 2: 3, 4: 4, 5: 5}
>>> events = [COAdsEvent, CODesEvent, OAdsEvent,
>>>           ODesEvent, CODiffEvent,
>>>           ODiffEvent, COOxEvent]

Here event 0 has a reverse event 1, 2 has 3, 4 and 5 are their own inverses because they are diffusion, and 6 is left out because it is assumed irreversible.
The numbering of events is determined by the order in the list :code:`events` defined here.

Define and run simulation
-----------------------------

Now the simulation object :class:`NeighborKMC.user_kmc.NeighborKMC` can be defined and the simulation performed:

>>> # Instantiate simulator object.
>>> sim = NeighborKMC(system=p, tend=tend,
>>>                   parameters=parameters,
>>>                   events=events,
>>>                   rev_events=reverse_events)
>>> result = sim.run_kmc()
>>> print("Simulation end time reached ! ! !")




.. _analyzecoox:

Analyze results
----------------------------
The results are analyzed by reading in the :ref:`code output <output>`. For example, if we need to calculate the CO and O coverage as a function of time for the entire system:

>>> import numpy as np
>>> time = np.loadtxt("time.txt")
>>> covs = np.loadtxt("coverages.txt")
>>> Nsites = float(len(covs[0]))
>>> cov_CO = [sum([1 for val in covs[i] if val == 1]) / Nsites for i in range(len(covs))]
>>> cov_O = [sum([1 for val in covs[i] if val == 2]) / Nsites for i in range(len(covs))]
>>> cov_free = [sum([1 for val in covs[i] if val == 0]) / Nsites for i in range(len(covs))]

If we need to analyze it for each site-type, the site-types need to be read. For corners, this may look like:

>>> stypes = np.loadtxt("stypes.txt")
>>> Ncnr = float(len([s for s in stypes if s == 0])) # Number of corner sites
>>> cov_CO_corners = [sum([1 for val in covs[i] if val == 1]) / Ncnr for i in range(len(covs))]

Typically, a turnover frequency is also relevant to calculate:

>>> Nevents = 7 # How many types of events are there.
>>> sid_ev = np.loadtxt("sid_ev.txt").reshape(-1,stypes.shape[0],Nevents)
>>> sid_ev_other = np.loadtxt("sid_ev.txt").reshape(-1,stypes.shape[0],Nevents)
>>> TOF = sum([sid_ev[-1][-1]+sid_ev_other[-1][-1]) / (Nsites*time[-1]) # How many CO+O->CO2 has fired per time and site. 

Often it can be useful to discard points out of steady-state by selecting only part of :code:`sid_ev`.
To draw statistically sound conclusions, it is recommended that multiple identically prepared simulations are performed.
See the tutorial on :ref:`Parallel simulations <parallel>`.





