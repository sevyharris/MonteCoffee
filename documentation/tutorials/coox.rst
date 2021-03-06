.. _coox:
.. index:: CO oxidation

CO oxidation on Pt(111)
*************************************

We present here the CO oxidation on Pt(111) which is published along with the oxidation on Pt nanoparticles in: `M. Jørgensen and H. Grönbeck, ACS Catal., 7, 5054-5061 (2017) <https://pubs.acs.org/doi/10.1021/acscatal.7b01194>`_ . 

For the CO oxidation the following five chemical reactions have to be considered:

.. math:: 
   
   O_2 + 2^* & \longleftrightarrow  2O^* \\
   CO + * & \longleftrightarrow  CO* \\
   CO^* + O^* & \longrightarrow  CO_2 \\
   CO^* + * & \longleftrightarrow  * + CO^* \\
   O^* + * & \longleftrightarrow  * + O^*

   
which are dissociative adsorption of oxygen, adsorption of CO, reaction of adsorbed O and CO to CO\ :sub:`2` and CO and O diffusion. In the model it is assumed that the reaction to CO\ :sub:`2` is not reversible and the reaction product is directly desorbing. 

In this example, species 0 denotes empty sites, 1 is CO and 2 refers to O. All energies were obtained using density functional theory and are given in the paper mentioned above. 

Define constants and parameters
-------------------------------

Any global constants can be defined in `user_constants.py <../api/NeighborKMC.examples.surface.html#module-NeighborKMC.examples.Pt_111_COOx.user_constants>`_. In this example this file stores physical constants but also the calculated vibrations of the various adsorbed species and of the CO\ :sub:`2` formation transition state and their physical properties.
 
Next, the simulation parameters must be set and stored in a dictionary:

.. code-block:: python

     T = 800.  # Temperature (K)
     pCO = 2E3  # CO pressure (Pa)
     pO2 = 1E3  # O2 pressure (Pa)
     tend = 1E-3  # End time of simulation (s)

Here the Temperature is set to 800 K, and the pressure ratio of CO to O\ :sub:`2` to 2:1. The pressure is given in Pa. Also the end time is defined, but because we are interested mainly in steady-state properties, the simulation is run until steady-state is reached and usually stopped before the given end time reached. 

Define sites and system
----------------------------

One site is defined for each surface atom using an :class:`Ase.Atoms` object. We are using the :class:`ase.build` class to construct our (111) surface. It consists of only one layer, because its only purpose is to simplify the site-connectivity set-up. Therefore, the used lattice-parameter is also not related to any results of density functional theory calculations or experiments. 
Thus starting with an empty list of sites we construct a :math:`10x10` surface consisting of 100 atoms:

.. code-block:: python

     from ase.build import fcc111
     from user_sites import Site

     sites = []
     a = 4.00  # Lattice Parameter (not related to DFT!)
     Ncutoff = a / np.sqrt(2.) + 0.05  # Nearest neighbor cutoff
     atoms = fcc111("Pt", a = a, size = (10,10,1))
     atoms.write('surface.traj')

We also have written out the prepared surface as `.traj` file to see if it looks as intended. As next step we have to add a site for each surface atom. Although we know that the oxygen preferably adsorbs on fcc sites and the CO on top sites, we assume here all sites to be equal (coarse-graining) but use the corresponding energies of the preferred sites.

.. code-block:: python

     # Create a site for each surface-atom:
     for i in range(len(atoms)): 
         sites.append(Site(stype=0,
                           covered=0, ind=i))

Here, the block :code:`ind=i` stores the index of the atom in the :class:`ASE.Atoms` object on the :class:`NeighborKMC.user_sites.Site` object, which can be later on used to write out `.traj`-files during the simulation. 

Now the :class:`NeighborKMC.user_system.System` object can be defined from the collection of sites:

.. code-block:: python

     from user_system import System
     p = System(atoms=atoms, # store ASE.Atoms as well
                sites=sites)


In addition to the overall system, we also need to define neighbor lists. It is simplest to define this according to the nearest neighbor distances. Thus we call in the main definition file of our kMC simulation:

.. code-block:: python

     p.set_neighbors(Ncutoff,pbc = True)

which sets the global neighbor list based on distances for our System. Because we are using a surface, periodic boundary conditions are desirable. The actual function is then defined in :class:`NeighborKMC.user_sites.Site` and looks as follows:

.. code-block:: python

    def set_neighbors(self, Ncutoff, pbc=False):

        if self.atoms is None:
            raise Warning("Tried to set neighbor-distances in user_system.set_neighbors() with self.atom = None")

        for i, s in enumerate(self.sites):  #For all sites
            for j, sother in enumerate(self.sites): #Check all the other sites 
                dcur = self.atoms.get_distance(s.ind, sother.ind, mic=pbc) # use ase function get_distance
                if dcur < Ncutoff and j != i:
                    s.neighbors.append(j)        #add site j to neighbor list of site i
        if len(self.neighbors) == 0:             #check if neighbors exists -- otherwise the site will not interact with each other
            self.neighbors = [s.neighbors for s in self.sites]
            self.verify_nlist()


 
Define reaction energies and entropies
--------------------------------------------
In this step, the reaction energies, or methods to calculate these, are defined in `user_energy.py <../api/NeighborKMC.examples.surface.html#module-NeighborKMC.examples.Pt_111_COOx.user_energy>`_ and the entropies in `user_entropy.py <../api/NeighborKMC.examples.surface.html#module-NeighborKMC.examples.Pt_111_COOx.user_entropy>`_. That makes it simple to do all the book-keeping accordingly. 

We used the from density functional theory obtained energies for the adsorption of oxygen on the fcc site and of CO on the top site and diffusion constants for oxygen to diffuse from fcc to fcc site and for CO for from top to top site. 
In this example, the adsorption energies are defined to be positiv in contradiction to the more generally used negative notation in theoretical papers. 

.. code-block:: python

     EadsCO = 1.36 
     EadsO = 0.97 

     EdiffCO = 0.08 # CO diffusion barrier
     EdiffO = 0.58 # O diffusion barrier

From the adsorption energies the activation energy for the CO\ :sub:`2` formation is calculated from the BEP relations according to: 

.. code-block:: python

    def get_Ea(ECO, EO):
        ETS = 0.824 * (-EO -ECO) + 0.168 + 0.47238  # How much larger is the energy of CO and O wrt Pt(111)
        Ea = ETS + ECO + EO  # Translate the barriers relative to Pt(111)
        return Ea

The reason why not one single activation energy is used, are the repulsive adsorbate-adsorbate interactions which depend locally on the coverage of the neighbors of the reaction site:

.. code-block:: python

     def get_repulsion(cov_self, cov_NN, stype):
         repulsion = 0.
         ECOCO = 0.19  #  How CO affects CO
         EOO = 0.32  # How O affects O 
         ECOO = 0.3  # How CO affects O
         EOCO = 0.3  # How O affects CO
         HInttwo = [[0., 0., 0.], [0., ECOCO, EOCO],
                   [0., ECOO, EOO]]  # Two body interaction Hamiltonian 3x3 because 0 = empty.
         for j in cov_NN:  # For each covered Neighbor, give a repulsion:
             repulsion += HInttwo[cov_self][j]

         return repulsion

The in `user_entropy.py <../api/NeighborKMC.examples.surface.html#module-NeighborKMC.examples.Pt_111_COOx.user_entropy>`_ stored entropies are included in the calculation of the rates as will be shown in the next section. Generally we define here the entropy for gas-phase CO and oxygen, as well as a method to calculate harmonic adsorbate entropy. The definitions of this functions can be looked up in this file directly. 

.. _defeventsquick:

Define events
--------------
Here event-types are defined, which are stored in `user_events.py <../api/NeighborKMC.examples.surface.html#module-NeighborKMC.examples.Pt_111_COOx.user_events>`_.
For each possible type of event, a class is derived from :class:`NeighborKMC.base.events.EventBase`.
For example take the :math:`\mathrm{CO_2}` formation. This event is defined in `user_events.py <../api/NeighborKMC.examples.surface.html#module-NeighborKMC.examples.Pt_111_COOx.user_events>`_ as follows.

First we import the necessary functions, classes, and constants:

.. code-block:: python

     from base.events import EventBase
     from user_entropy import get_Zvib, get_Z_CO, get_Z_O2
     from user_constants import mCO, mO2, Asite, modes_COads, modes_Oads,\
             modes_TS_COOx, modes_COgas, modes_O2gas, kB, eV2J, s0CO, s0O, h,
     from user_energy import EadsCO, EadsO, get_Ea, get_repulsion, EdiffCO, EdiffO

Now we derive a class to contain the event:

.. code-block:: python

     class COOxEvent(EventBase):
         def __init__(self, params):
             Zads = get_Zvib(params["T"], modes_COads) * get_Zvib(params["T"], modes_Oads)
             Zts = get_Zvib(params["T"], modes_TS_COOx)
             self.Zratio = Zts / Zads 
             EventBase.__init__(self, params)

The constructor :code:`__init__(self,params)`  attaches relevant parameters to the object, and :code:`self.Zratio` is the ratio
between the partition functions in the initial state and transition state, used to calculate the rate constant in transition state theory. We need a function `possible(self,system, site, other_site)`
that returns True if the event is possible on the current site-pair:

.. code-block:: python

         def possible(self, system, site, other_site):
             if (system.sites[site].covered == 1 and system.sites[other_site].covered == 2) or \
                   (system.sites[site].covered == 2 and system.sites[other_site].covered == 1):
                 return True
             else:
                 return False

Thus, for the event to be possible, the site needs to be covered by 1 (CO) and the neighbor site by 2 (O) or the other way round. That is originated by the use of single neighbor site pairs. Thus a pair with the indexes (10,11) would be the same as (11,10) in the code to avoid double counting in the time list.

Next we need to define a function :code:`get_rate(self, system, i_site, other_site)` that returns the rate constant:

.. code-block:: python

        def get_rate(self, system, i_site, other_site):
            ECO = EadsCO
            EO = EadsO
            if system.sites[site].covered == 1:
                Ncovs_CO = [system.sites[n].covered for n in system.neighbors[site] ]
                Ncovs_O = [system.sites[n].covered for n in system.neighbors[other_site]]
            else:
                Ncovs_CO = [system.sites[n].covered for n in system.neighbors[other_site] ]
                Ncovs_O = [system.sites[n].covered for n in system.neighbors[site]]
            ECO -= get_repulsion(1, Ncovs_CO, 0)
            EO -= get_repulsion(2, Ncovs_O, 0)
            Ea = max(0., get_Ea(ECO, EO)) # No negative energy barriers
            return self.alpha * self.Zratio * np.exp(-Ea/(kB * self.params['T'])) * kB * self.params['T'] / h

The adsorption energies are fixed, but weaken in the case of any repulsions. 
A call is made to :code:`get_Ea(ECO, EO)` to obtain the reaction energy barrier.
The rate constant is multiplied with :code:`self.alpha`, because :code:`self.alpha` is the slowing-down factor that is adjusted dynamically for each event during simulation.

Also the event requires a method :code:`do_event(self,system, site, other_site)` to perform modifications to the site-occupations when fired:

.. code-block:: python

        def do_event(self, system, site, other_site):
            system.sites[site].covered = 0
            system.sites[other_site].covered = 0

In this case, the two sites containing CO and O are simply emptied. At the end we define the method :code:`get_involve_other(self)` to define if the neighboring sites are actually involved in the event:

.. code-block:: python

        def get_involve_other(self):
            return True

After giving this example for one event, the other events can be defined similarly. How to define single site and dissociative adsorption was shown before. Only the rates have to be adjusted according to transition state theory. Having all events for each type of reaction defined in this order:

    - (0) :class:`NeighborKMC.user_events.COAdsEvent` for CO adsorption.
    - (1) :class:`NeighborKMC.user_events.CODesEvent` for CO desorption.
    - (2) :class:`NeighborKMC.user_events.OAdsEvent` for O2 dissociative adsorption.
    - (3) :class:`NeighborKMC.user_events.ODesEvent` for O2 desorption.
    - (4) :class:`NeighborKMC.user_events.CODiffEvent` for CO diffusion.
    - (5) :class:`NeighborKMC.user_events.ODiffEvent` for O diffusion.
    - (6) :class:`NeighborKMC.user_events.COOxEvent` for CO+O -> CO2.

we can now store the event-class references in a list for the NeighborKMC simulations and define accordingly a list of reverse events. In this example the CO (O) adsorption and desorption are reverse to each other. In addition are the diffusion processes reverse to themselves. The CO oxidation itself is not reversible. Thus it isn't defined in this list.

.. code-block:: python

     events = [COAdsEvent, CODesEvent, OAdsEvent,
               ODesEvent, CODiffEvent,
               ODiffEvent, COOxEvent]
     reverse_events = {0: 1, 2: 3, 4: 4, 5: 5}

The numbering of events is determined by the order in the list :code:`events` defined here.

Define and run simulation
-----------------------------

Now the simulation object :class:`NeighborKMC.user_kmc.NeighborKMC` can be defined and the simulation performed:

.. code-block:: python

     # Instantiate simulator object.
     sim = NeighborKMC(system=p, tend=tend,
                       parameters=parameters,
                       events=events,
                       rev_events=reverse_events)
     sim.run_kmc()
     print("Simulation end time reached ! ! !")




.. _analyzecoox:

Analyze results
----------------------------
The results are analyzed by reading in the :ref:`code output <output>`. For example, if we need to calculate the CO and O coverage as a function of time for the entire system:

.. code-block:: python

     import numpy as np
     time = np.loadtxt("time.txt")
     covs = np.loadtxt("coverages.txt")
     Nsites = float(len(covs[0]))
     cov_CO = [sum([1 for val in covs[i] if val == 1]) / Nsites for i in range(len(covs))]
     cov_O = [sum([1 for val in covs[i] if val == 2]) / Nsites for i in range(len(covs))]
     cov_free = [sum([1 for val in covs[i] if val == 0]) / Nsites for i in range(len(covs))]

Typically, a turnover frequency is also relevant to calculate:

.. code-block:: python

     evs_exec = np.loadtxt("evs_exec.txt")
     TOF = evs_exec[-1] / (Nsites*time[-1]) # How many CO+O->CO2 has fired per time and site.

Often it can be useful to discard points out of steady-state. The detailed evolution of the time with the sites can be found in `detail_site_event_evol.hdf5`.
To draw statistically sound conclusions, it is recommended that multiple identically prepared simulations are performed and in this case at least 100 CO oxidation events are performed in the steady state. 



